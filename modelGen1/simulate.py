#!/usr/bin/env python3
"""OpenPhase 多晶 PVD：以 vapor->grain 成对驱动力耦合定向沉积通量。"""
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

# 这些常量必须与 ProjectInput.opi 保持一致。当前 OpenPhase Python binding
# 没有完整暴露 GridParameters，因此输出 VTK 时在脚本中显式记录网格尺寸。
NX, NY, DX, DT = 192, 128, 2.0e-9, 1.0e-8
VAPOR_PHASE, SOLID_PHASE = 0, 1
GRAIN_LABEL_THRESHOLD = 0.50
SEED_RADIUS_GRID = 8.5
SEED_CENTERS_X = tuple(float(x) for x in range(8, NX, 16))


def import_openphase():
    """导入 OpenPhase；若沉积桥接模块不存在，则先调用 Makefile 编译它。"""
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
        import openphase_deposition as coupling
    except ImportError:
        print("Building local OpenPhase deposition coupling ...", flush=True)
        subprocess.run(["make"], cwd=CASE_DIR, check=True)
        import openphase_deposition as coupling
    return op, coupling


def parse_args() -> argparse.Namespace:
    """读取命令行参数；这些参数只改变外部沉积条件和输出控制。"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--steps", type=int, default=20_000)
    parser.add_argument("--output-every", type=int, default=500)
    parser.add_argument("--output", type=Path, default=CASE_DIR / "output")
    parser.add_argument(
        "--flux", type=float, default=1.0,
        help="入射通量幅值；第一阶段使用归一化单位（默认 1）",
    )
    parser.add_argument(
        "--incidence-angle", type=float, default=0.0,
        help="相对基底法线的入射角（度）；正值表示通量具有 +x 分量",
    )
    parser.add_argument(
        "--flux-to-driving", type=float, default=2.0e7,
        help="通量到峰值体积驱动力的换算系数，单位 Pa/flux-unit",
    )
    parser.add_argument(
        "--partition-power", type=float, default=4.0,
        help="三重点处按 phi_i**p 分配沉积驱动力的指数",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    """在创建 OpenPhase 对象前检查会导致数值或几何错误的参数。"""
    if args.steps < 0:
        raise SystemExit("--steps must be non-negative")
    if args.output_every <= 0:
        raise SystemExit("--output-every must be positive")
    if args.flux < 0.0:
        raise SystemExit("--flux must be non-negative")
    if args.flux_to_driving < 0.0:
        raise SystemExit("--flux-to-driving must be non-negative")
    if args.partition_power <= 0.0:
        raise SystemExit("--partition-power must be positive")
    if abs(args.incidence_angle) >= 90.0:
        raise SystemExit("--incidence-angle must lie strictly between -90 and 90 degrees")


def write_vtk(
    path: Path,
    vapor: np.ndarray,
    solid: np.ndarray,
    deposition_drive: np.ndarray,
    grain_id: np.ndarray,
) -> None:
    """写出 ParaView 可读的 Legacy VTK；所有场均按 (y, x) 展平。"""
    fields = {
        "vapor_fraction": vapor,
        "solid_fraction": solid,
        "deposition_drive_Pa": deposition_drive,
    }
    for name, values in fields.items():
        if values.shape != (NY, NX) or not np.isfinite(values).all():
            raise ValueError(f"{name} must have shape {(NY, NX)} and finite values")
    if grain_id.shape != (NY, NX):
        raise ValueError(f"grain_id must have shape {(NY, NX)}")

    with path.open("w", encoding="ascii") as stream:
        stream.write("# vtk DataFile Version 3.0\n")
        stream.write("OpenPhase polycrystalline PVD\nASCII\n")
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

    # 这些对象分别负责输入参数、时间步、序参量、双障碍势、界面性质、
    # 成对驱动力和边界条件。真正的 MPF 离散更新由 OpenPhase 完成。
    settings = op.Settings(str(input_path))
    runtime = op.RunTimeControl(settings, str(input_path))
    phase = op.PhaseField(settings, str(input_path))
    obstacle = op.DoubleObstacle(settings, str(input_path))
    interface = op.InterfaceProperties(settings, str(input_path))
    driving = op.DrivingForce(settings, str(input_path))
    boundary = op.BoundaryConditions(settings, str(input_path))

    # OpenPhase 的 phase index 0 是一个气相序参量；phase index 1 对应
    # 热力学上的 Film 相，但每次 Sphere 初始化会创建一个独立的 grain
    # 序参量。因此这里得到“1 个 vapor + 多个 grain”的多相场系统。
    vapor_grain = op.Initializations.Single(phase, VAPOR_PHASE, boundary)
    phase.FieldsProperties[vapor_grain].State = op.AggregateStates.Gas
    grains: list[int] = []
    for center_x in SEED_CENTERS_X:
        # 半径以网格点计，圆心放在 y=0，OpenPhase 的 sectional sphere
        # 形成贴附基底的半圆晶种；晶种沿周期 x 方向相切排列。
        grain = op.Initializations.Sphere(
            phase, SOLID_PHASE, SEED_RADIUS_GRID,
            center_x, 0.0, 0.0, boundary, True,
        )
        phase.FieldsProperties[grain].State = op.AggregateStates.Solid
        grains.append(int(grain))

    # 入射角相对基底法线定义：0° 为竖直向下，正角度增加 +x 分量。
    # 向量方向从沉积源指向薄膜，因此 flux_y 为负值。
    angle = math.radians(args.incidence_angle)
    flux_x = args.flux * math.sin(angle)
    flux_y = -args.flux * math.cos(angle)

    # 将本次运行的控制参数保存下来，避免只凭 VTK 文件无法复现实验。
    with (args.output / "run_parameters.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(("parameter", "value"))
        writer.writerows((
            ("steps", args.steps),
            ("dt_s", DT),
            ("dx_m", DX),
            ("flux_normalized", args.flux),
            ("incidence_angle_deg", args.incidence_angle),
            ("flux_x", flux_x),
            ("flux_y", flux_y),
            ("flux_to_driving_Pa_per_flux_unit", args.flux_to_driving),
            ("partition_power", args.partition_power),
            ("number_of_grains", len(grains)),
        ))

    # 计算初始界面上的沉积驱动力场。它只用于 step=0 输出；每个时间步
    # 会在加入当前 MPF 曲率驱动力后重新计算。
    deposition = coupling.evaluate_deposition(
        phase, VAPOR_PHASE, flux_x, flux_y,
        args.flux_to_driving, args.partition_power,
    )

    with (args.output / "history.csv").open(
        "w", newline="", encoding="utf-8"
    ) as history_stream:
        history = csv.writer(history_stream)
        history.writerow((
            "step", "time_s", "mean_solid_fraction",
            "equivalent_film_thickness_m", "maximum_height_m",
            "maximum_partition_error", "maximum_deposition_drive_Pa",
        ))

        def output(step: int) -> None:
            # 聚合所有 solid grain 得到固相分数；dominant_grain 只用于
            # 可视化标签，不参与相场演化。弥散界面处标签设为 -1。
            vapor = coupling.phase_fraction(phase, VAPOR_PHASE)
            solid = coupling.phase_fraction(phase, SOLID_PHASE)
            grain_id = coupling.dominant_grain(phase, SOLID_PHASE)
            grain_id[solid < GRAIN_LABEL_THRESHOLD] = -1
            write_vtk(
                args.output / f"fields_{step:06d}.vtk",
                vapor, solid, deposition, grain_id,
            )
            # 最高高度是网格阈值统计；等效高度用固相体积/基底宽度，
            # 对弥散界面更平滑，适合比较成膜速率。
            occupied_rows = np.where(np.any(solid >= GRAIN_LABEL_THRESHOLD, axis=1))[0]
            max_height = (occupied_rows.max() + 1) * DX if occupied_rows.size else 0.0
            equivalent_height = float(solid.sum()) * DX / NX
            partition_error = float(coupling.maximum_partition_error(phase))
            max_drive = float(np.max(deposition))
            history.writerow((
                step, f"{step * DT:.9e}", f"{solid.mean():.9e}",
                f"{equivalent_height:.9e}", f"{max_height:.9e}",
                f"{partition_error:.9e}", f"{max_drive:.9e}",
            ))
            history_stream.flush()
            print(
                f"step={step:6d} t={step * DT:.3e} s "
                f"eq_height={equivalent_height * 1e9:7.2f} nm "
                f"max_height={max_height * 1e9:7.2f} nm "
                f"sum(phi)-1={partition_error:.2e}",
                flush=True,
            )

        output(0)
        for step in range(1, args.steps + 1):
            # 1) 清空上一步的成对驱动力。
            driving.Clear()
            # 2) 根据当前界面计算 OpenPhase 双障碍势/曲率贡献。
            interface.Set(phase, boundary)
            obstacle.CalculatePhaseFieldIncrements(phase, interface, driving)
            # 3) 计算局部法向投影 q，并把 vapor->grain 沉积项写入
            #    DrivingForce。多晶三重点按 phi_i**p 分配。
            deposition = coupling.add_deposition_driving_force(
                driving, phase, VAPOR_PHASE, SOLID_PHASE,
                flux_x, flux_y, args.flux_to_driving,
                args.partition_power,
            )
            # 4) 初始化/平均 pair weight 后合并所有驱动力，得到 d(phi)/dt。
            driving.Average(phase, boundary)
            driving.MergePhaseFieldIncrements(phase, interface)
            # 5) 归一化增量并推进相场；OpenPhase 的反对称 pair 更新
            #    保持 vapor + 所有 solid grain 的分数和为 1。
            phase.NormalizeIncrements(boundary, DT)
            phase.MergeIncrements(boundary, DT)
            runtime.IncrementTimeStep()
            if step % args.output_every == 0 or step == args.steps:
                output(step)

    print(f"Finished. Results are in {args.output.resolve()}")


if __name__ == "__main__":
    main()
