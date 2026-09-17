#!/usr/bin/env python3
"""OpenPhase 多晶 PVD：耦合 vapor density 对流、扩散和沉积消耗。"""
from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

CASE_DIR = Path(__file__).resolve().parent
OPENPHASE_DIR = Path("/home/majie/softwares/OpenPhase")
OPENPHASE_BINDINGS = OPENPHASE_DIR / "pythonbindings"
DEFAULT_INPUT = CASE_DIR / "ProjectInput.opi"

# 必须与 ProjectInput.opi 保持一致；当前 binding 未完整暴露网格参数。
NX, NY, DX, DT = 192, 128, 2.0e-9, 1.0e-8
VAPOR_PHASE, SOLID_PHASE = 0, 1
GRAIN_LABEL_THRESHOLD = 0.50
SEED_RADIUS_GRID = 8.5
SEED_CENTERS_X = tuple(float(x) for x in range(8, NX, 16))


def import_openphase():
    """导入 OpenPhase；缺少本地 vapor 耦合模块时自动编译。"""
    sys.path.insert(0, str(OPENPHASE_BINDINGS))
    sys.path.insert(0, str(CASE_DIR))
    try:
        import OpenPhase as op
    except ImportError as exc:
        raise SystemExit(
            f"Cannot import OpenPhase from {OPENPHASE_BINDINGS}. Run with "
            "/home/majie/trainFolder/openphaseFolder/op_env/bin/python."
        ) from exc
    try:
        import openphase_vapor as coupling
    except ImportError:
        print("Building local OpenPhase vapor coupling ...", flush=True)
        subprocess.run(["make"], cwd=CASE_DIR, check=True)
        import openphase_vapor as coupling
    return op, coupling


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--steps", type=int, default=20_000)
    parser.add_argument("--output-every", type=int, default=500)
    parser.add_argument("--output", type=Path, default=CASE_DIR / "output")
    parser.add_argument("--flux", type=float, default=1.0,
                        help="归一化入射粒子通量")
    parser.add_argument("--incidence-angle", type=float, default=0.0,
                        help="相对基底法线的入射角（度），正值朝 +x 倾斜")
    parser.add_argument("--flux-to-driving", type=float, default=2.0e7,
                        help="rho*flux 到体积驱动力的换算系数 [Pa/flux-unit]")
    parser.add_argument("--partition-power", type=float, default=4.0,
                        help="三重点沉积驱动力按 phi_i**p 分配的指数")
    parser.add_argument("--vapor-feed-density", type=float, default=1.0,
                        help="顶部 Dirichlet vapor density")
    parser.add_argument("--vapor-diffusivity", type=float, default=2.0e-11,
                        help="vapor diffusivity [m^2/s]")
    parser.add_argument("--transport-speed", type=float, default=2.0e-2,
                        help="开启对流时的速度幅值 [m/s]")
    parser.add_argument(
        "--convection", action=argparse.BooleanOptionalAction, default=True,
        help="开启/关闭 vapor density 的对流项（可用 --no-convection 关闭）",
    )
    parser.add_argument("--consumption-factor", type=float, default=1.0,
                        help="单位固相增量消耗的归一化 vapor density")
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.steps < 0 or args.output_every <= 0:
        raise SystemExit("--steps must be non-negative and --output-every positive")
    nonnegative = {
        "--flux": args.flux,
        "--flux-to-driving": args.flux_to_driving,
        "--vapor-feed-density": args.vapor_feed_density,
        "--vapor-diffusivity": args.vapor_diffusivity,
        "--transport-speed": args.transport_speed,
        "--consumption-factor": args.consumption_factor,
    }
    for name, value in nonnegative.items():
        if value < 0.0:
            raise SystemExit(f"{name} must be non-negative")
    if args.partition_power <= 0.0:
        raise SystemExit("--partition-power must be positive")
    if abs(args.incidence_angle) >= 90.0:
        raise SystemExit("--incidence-angle must lie strictly between -90 and 90 degrees")


def transport_substeps(diffusivity: float, vx: float, vy: float) -> int:
    """按显式扩散和一阶迎风 CFL 数自动选择 transport 子步数。"""
    diffusion_number = 4.0 * diffusivity * DT / DX**2
    advection_number = (abs(vx) + abs(vy)) * DT / DX
    return max(1, math.ceil((diffusion_number + advection_number) / 0.8))


def advance_vapor_density(
    density: np.ndarray,
    vapor_fraction: np.ndarray,
    diffusivity: float,
    vx: float,
    vy: float,
    feed_density: float,
) -> tuple[np.ndarray, int]:
    """用守恒面通量推进带气相掩膜的二维对流—扩散方程。

    x 方向周期；顶部固定供料；底部无通量。对流采用一阶迎风格式，
    扩散采用中心差分。面通量乘相邻单元气相分数的较小值，使 transport
    不穿过致密固相。返回更新后的密度和实际使用的子步数。
    """
    substeps = transport_substeps(diffusivity, vx, vy)
    sub_dt = DT / substeps
    rho = density.copy()
    access = np.clip(vapor_fraction, 0.0, 1.0)

    # 气相掩膜及其四个面上的可通行比例。
    access_right = np.minimum(access, np.roll(access, -1, axis=1))
    access_left = np.roll(access_right, 1, axis=1)
    access_top = np.empty_like(access)
    access_top[:-1] = np.minimum(access[:-1], access[1:])
    access_top[-1] = access[-1]
    access_bottom = np.zeros_like(access)
    access_bottom[1:] = access_top[:-1]

    for _ in range(substeps):
        right = np.roll(rho, -1, axis=1)
        left = np.roll(rho, 1, axis=1)
        top = np.empty_like(rho)
        top[:-1] = rho[1:]
        top[-1] = feed_density
        bottom = np.empty_like(rho)
        bottom[1:] = rho[:-1]
        bottom[0] = rho[0]  # 基底为无通量边界

        diffusion = diffusivity * (
            access_right * (right - rho) + access_left * (left - rho)
            + access_top * (top - rho) + access_bottom * (bottom - rho)
        ) / DX**2

        # x 面的迎风对流通量；正号表示沿 +x 方向。
        flux_right = vx * (rho if vx >= 0.0 else right) * access_right
        flux_left = np.roll(flux_right, 1, axis=1)

        # y 顶面通量。vy<0 时上方单元是 donor，最顶部 donor 为供料边界。
        flux_top = vy * (rho if vy >= 0.0 else top) * access_top
        flux_bottom = np.zeros_like(rho)
        flux_bottom[1:] = flux_top[:-1]
        advection = -(
            flux_right - flux_left + flux_top - flux_bottom
        ) / DX

        rho += sub_dt * (diffusion + advection)
        rho = np.clip(rho, 0.0, feed_density)
        rho[vapor_fraction < 1.0e-6] = 0.0
        rho[-1, :] = feed_density
    return rho, substeps


def write_vtk(
    path: Path,
    vapor: np.ndarray,
    solid: np.ndarray,
    vapor_density: np.ndarray,
    deposition_drive: np.ndarray,
    grain_id: np.ndarray,
) -> None:
    """写出 ParaView 可读的 Legacy VTK。"""
    fields = {
        "vapor_fraction": vapor,
        "solid_fraction": solid,
        "vapor_density": vapor_density,
        "deposition_drive_Pa": deposition_drive,
    }
    for name, values in fields.items():
        if values.shape != (NY, NX) or not np.isfinite(values).all():
            raise ValueError(f"{name} must have shape {(NY, NX)} and finite values")
    with path.open("w", encoding="ascii") as stream:
        stream.write("# vtk DataFile Version 3.0\n")
        stream.write("OpenPhase PVD with vapor transport\nASCII\n")
        stream.write("DATASET STRUCTURED_POINTS\n")
        stream.write(f"DIMENSIONS {NX} {NY} 1\nORIGIN 0 0 0\n")
        stream.write(f"SPACING {DX:.12e} {DX:.12e} {DX:.12e}\n")
        stream.write(f"POINT_DATA {NX * NY}\n")
        for name, values in fields.items():
            stream.write(f"SCALARS {name} double 1\nLOOKUP_TABLE default\n")
            np.savetxt(stream, values.reshape(-1), fmt="%.10e")
        stream.write("SCALARS grain_id int 1\nLOOKUP_TABLE default\n")
        np.savetxt(stream, grain_id.reshape(-1), fmt="%d")


def main() -> None:
    args = parse_args()
    validate_args(args)
    op, coupling = import_openphase()
    input_path = args.input.resolve()
    args.output.mkdir(parents=True, exist_ok=True)

    settings = op.Settings(str(input_path))
    runtime = op.RunTimeControl(settings, str(input_path))
    phase = op.PhaseField(settings, str(input_path))
    obstacle = op.DoubleObstacle(settings, str(input_path))
    interface = op.InterfaceProperties(settings, str(input_path))
    driving = op.DrivingForce(settings, str(input_path))
    boundary = op.BoundaryConditions(settings, str(input_path))

    # 一个 vapor grain 填充计算域，多个 Film grain 以半圆晶种贴在 y=0。
    vapor_grain = op.Initializations.Single(phase, VAPOR_PHASE, boundary)
    phase.FieldsProperties[vapor_grain].State = op.AggregateStates.Gas
    grains: list[int] = []
    for center_x in SEED_CENTERS_X:
        grain = op.Initializations.Sphere(
            phase, SOLID_PHASE, SEED_RADIUS_GRID,
            center_x, 0.0, 0.0, boundary, True,
        )
        phase.FieldsProperties[grain].State = op.AggregateStates.Solid
        grains.append(int(grain))

    angle = math.radians(args.incidence_angle)
    flux_x = args.flux * math.sin(angle)
    flux_y = -args.flux * math.cos(angle)
    if args.convection:
        velocity_x = args.transport_speed * math.sin(angle)
        velocity_y = -args.transport_speed * math.cos(angle)
    else:
        velocity_x = velocity_y = 0.0

    vapor = coupling.phase_fraction(phase, VAPOR_PHASE)
    solid = coupling.phase_fraction(phase, SOLID_PHASE)
    # density 是气相中的归一化局部密度；致密固相内部初始为零。
    vapor_density = np.full((NY, NX), args.vapor_feed_density, dtype=float)
    vapor_density[vapor < 1.0e-6] = 0.0
    vapor_density[-1, :] = args.vapor_feed_density
    deposition = coupling.evaluate_deposition(
        phase, vapor_density, VAPOR_PHASE, flux_x, flux_y,
        args.flux_to_driving, args.partition_power,
    )
    substeps = transport_substeps(args.vapor_diffusivity, velocity_x, velocity_y)

    with (args.output / "run_parameters.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(("parameter", "value"))
        writer.writerows((
            ("steps", args.steps), ("dt_s", DT), ("dx_m", DX),
            ("flux_normalized", args.flux),
            ("incidence_angle_deg", args.incidence_angle),
            ("flux_to_driving_Pa_per_flux_unit", args.flux_to_driving),
            ("partition_power", args.partition_power),
            ("vapor_feed_density", args.vapor_feed_density),
            ("vapor_diffusivity_m2_per_s", args.vapor_diffusivity),
            ("convection", args.convection),
            ("transport_speed_m_per_s", args.transport_speed),
            ("velocity_x_m_per_s", velocity_x),
            ("velocity_y_m_per_s", velocity_y),
            ("transport_substeps", substeps),
            ("consumption_factor", args.consumption_factor),
            ("number_of_grains", len(grains)),
        ))

    cumulative_consumption = 0.0
    with (args.output / "history.csv").open(
        "w", newline="", encoding="utf-8"
    ) as history_stream:
        history = csv.writer(history_stream)
        history.writerow((
            "step", "time_s", "mean_solid_fraction",
            "equivalent_film_thickness_m", "maximum_height_m",
            "mean_vapor_density", "total_vapor_density",
            "cumulative_deposition_consumption",
            "maximum_partition_error", "maximum_deposition_drive_Pa",
            "transport_substeps",
        ))

        def output(step: int) -> None:
            grain_id = coupling.dominant_grain(phase, SOLID_PHASE)
            grain_id[solid < GRAIN_LABEL_THRESHOLD] = -1
            write_vtk(
                args.output / f"fields_{step:06d}.vtk",
                vapor, solid, vapor_density, deposition, grain_id,
            )
            occupied = np.where(np.any(solid >= GRAIN_LABEL_THRESHOLD, axis=1))[0]
            max_height = (occupied.max() + 1) * DX if occupied.size else 0.0
            equivalent_height = float(solid.sum()) * DX / NX
            vapor_cells = vapor > 0.5
            mean_density = (
                float(vapor_density[vapor_cells].mean())
                if np.any(vapor_cells) else 0.0
            )
            partition_error = float(coupling.maximum_partition_error(phase))
            history.writerow((
                step, f"{step * DT:.9e}", f"{solid.mean():.9e}",
                f"{equivalent_height:.9e}", f"{max_height:.9e}",
                f"{mean_density:.9e}", f"{vapor_density.sum():.9e}",
                f"{cumulative_consumption:.9e}", f"{partition_error:.9e}",
                f"{np.max(deposition):.9e}", substeps,
            ))
            history_stream.flush()
            print(
                f"step={step:6d} eq_height={equivalent_height * 1e9:7.2f} nm "
                f"rho_mean={mean_density:.4f} "
                f"sum(phi)-1={partition_error:.2e}", flush=True,
            )

        output(0)
        for step in range(1, args.steps + 1):
            old_solid = solid.copy()
            driving.Clear()
            interface.Set(phase, boundary)
            obstacle.CalculatePhaseFieldIncrements(phase, interface, driving)
            # 局部驱动力正比于 rho*|J|，并只作用于迎流气固界面。
            deposition = coupling.add_deposition_driving_force(
                driving, phase, vapor_density, VAPOR_PHASE, SOLID_PHASE,
                flux_x, flux_y, args.flux_to_driving, args.partition_power,
            )
            driving.Average(phase, boundary)
            driving.MergePhaseFieldIncrements(phase, interface)
            phase.NormalizeIncrements(boundary, DT)
            phase.MergeIncrements(boundary, DT)

            vapor = coupling.phase_fraction(phase, VAPOR_PHASE)
            solid = coupling.phase_fraction(phase, SOLID_PHASE)
            # 只有沉积引起的正固相增量消耗 vapor；曲率导致的溶解不返还。
            solid_increment = np.maximum(solid - old_solid, 0.0)
            local_consumption = args.consumption_factor * solid_increment
            cumulative_consumption += float(local_consumption.sum())
            vapor_density = np.clip(
                vapor_density - local_consumption,
                0.0, args.vapor_feed_density,
            )
            vapor_density, substeps = advance_vapor_density(
                vapor_density, vapor, args.vapor_diffusivity,
                velocity_x, velocity_y, args.vapor_feed_density,
            )
            runtime.IncrementTimeStep()
            if step % args.output_every == 0 or step == args.steps:
                output(step)

    print(f"Finished. Results are in {args.output.resolve()}")


if __name__ == "__main__":
    main()
