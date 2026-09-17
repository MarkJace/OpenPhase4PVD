# OpenPhase 多晶 PVD 相场模型：第一阶段技术说明

## 1. 本阶段的目标与边界

本目录实现薄膜沉积模型的第一步：只求解 OpenPhase 的 Steinbach 多相场，并把外部 deposition flux 转换为气相—固相界面的附加驱动力。模型包含多个独立固相晶粒，能描述气相被消耗、晶粒向上生长、晶界形成和晶粒竞争。

本阶段刻意不求解气相密度的对流—扩散，不包含 SiMTRA 粒子输运、表面扩散、组分、弹性、再成核和晶体各向异性。入射通量在全计算域中是给定的常向量。因此该模型是验证“多相场 + 沉积项”耦合是否正确的最小基线，不是已经标定的具体材料工艺模型。

参考文献为本目录上一级的 `ref_paper.pdf`：J. M. Monti et al., *Acta Materialia* 245 (2023) 118581。`path.md` 中关于把 vapor 作为一个相，并将沉积写成 pairwise driving force 的方案也被采用。

## 2. 从论文模型到 Steinbach/OpenPhase 模型

Monti 等人的多序参量模型写为

$$
\frac{\partial\phi_i}{\partial t}
=\frac{1}{\psi_s}\nabla\!\cdot\!\left[
\psi_s M(\phi_i)\nabla\frac{\delta\mathcal F}{\delta\phi_i}
\right]
+S(\boldsymbol\phi,\mathbf v)
\frac{\phi_i}{\sum_j\phi_j},
$$

其中论文式 (7) 的界面源项为

$$
S(\boldsymbol\phi,\mathbf v)=
\max\left\{-\left(\sum_i\phi_i^4\right)
(\rho\mathbf v)\cdot\nabla\psi_v,0\right\}.
$$

这个模型中的各晶粒序参量不是 OpenPhase 意义下局部和恒为 1 的一组相分数。Steinbach 多相场则要求

$$
\phi_v+\sum_{i=1}^{N}\phi_i=1.
$$

若直接给每个固相方程增加正源项，将破坏该约束。这里将沉积解释为局部的相分数转移

$$
\text{vapor}\rightarrow\text{grain }i,
$$

并把它写入 OpenPhase 的反对称成对驱动力容器。气相减少量和固相增加量天然相等，`NormalizeIncrements` 再处理有限时间步造成的数值越界，所以 partition of unity 得到保持。

## 3. 沉积驱动力

聚合气相分数为 $\phi_v$，从薄膜指向气相的局部法向定义为

$$
\mathbf n=\frac{\nabla\phi_v}{|\nabla\phi_v|}.
$$

输入的 $\mathbf J$ 是粒子从靶源指向薄膜的通量向量。局部朝向因子为

$$
q(\mathbf x)=\max\left(-\widehat{\mathbf J}\cdot\mathbf n,0\right).
$$

因此面向入射束的界面取 $q=1$，背向入射束的界面取 $q=0$。峰值沉积驱动力为

$$
\Delta G_0=\Lambda |\mathbf J|,
$$

其中 $\Lambda$ 是 `--flux-to-driving`。在气相和多个固相晶粒同时存在的弥散三重点，驱动力按

$$
w_i=\frac{\phi_i^p}{\sum_{j\in\mathrm{solid}}\phi_j^p},\qquad
\Delta G_{v i}^{\mathrm{dep}}=\Delta G_0q w_i
$$

分配，默认 $p=4$，与论文用四次幂抑制晶界处伪沉积的思想一致。在普通气固界面只有一个固相晶粒，故 $w_i=1$。OpenPhase 随后用其双障碍势 MPF 动力学、局部迁移率和 $\sqrt{\phi_v\phi_i}$ 插值把体积驱动力转换为界面运动。

需要区分两种“阴影”：当前 $q$ 能排除背向束流的局部表面，复现论文源项中的法向投影；它没有沿射线检查上游遮挡，所以深沟槽的全局 line-of-sight shadowing 要在后续加入射线追踪、可见度场或显式 vapor transport 后才完整。

## 4. 数值系统

计算为二维 $x$-$y$ 截面，$x$ 是基底方向，$+y$ 是薄膜生长方向。

| 项目 | 默认值 | 说明 |
|---|---:|---|
| 网格 | $192\times128$ | 物理尺寸 $384\times256$ nm |
| 网格间距 | 2 nm | 第一阶段兼顾界面分辨率与速度 |
| 弥散界面宽度 | 5 个网格 | OpenPhase `IWidth` |
| 时间步 | $10^{-8}$ s | 保守的显式 MPF 步长 |
| 初始晶粒 | 12 个半圆晶种 | 沿 $y=0$ 等间距排列并跨周期边界连续 |
| 气固界面能 | 1.0 J/m² | 原型参数 |
| 晶界能 | 0.50 J/m² | 原型参数 |
| 气固迁移率 | $5\times10^{-12}$ m⁴/(J·s) | 控制沉积界面响应 |
| 晶界迁移率 | $5\times10^{-14}$ m⁴/(J·s) | 允许较慢的晶粒竞争 |
| 峰值驱动力 | $2\times10^7$ Pa | `flux=1` 时的默认值 |
| 分配指数 | 4 | 三重点沉积分配 |

$x$ 两侧为周期边界；顶部和底部为无通量边界。底部半圆只在 $y\geq0$ 的计算域内出现，相当于贴附在平坦基底上的晶种。本阶段没有单独的基底相，因而底边只是不可穿越的几何边界。

这些参数优先保证数值原型可运行。它们不能直接解释为某种材料的定量参数。尤其 `flux-to-driving` 把归一化通量映射为 J/m³（即 Pa），应使用实验膜厚—时间曲线或原子体积、黏附系数和界面动力学进行标定。

## 5. 文件与软件结构

- `simulate.py`：建立多晶初态、推进 OpenPhase、输出 VTK 和历史量。
- `ProjectInput.opi`：OpenPhase 网格、相、界面能、迁移率、驱动力限制和边界条件。
- `openphase_deposition.cpp`：很小的 pybind11 桥接层；读取 OpenPhase 相场、计算 $q$ 和 $w_i$，并向 `DrivingForce::Force` 写入反对称气固驱动力。
- `Makefile`：编译桥接层并链接 `/home/majie/softwares/OpenPhase/lib/libOpenPhase.so`。
- `output/fields_*.vtk`：`vapor_fraction`、`solid_fraction`、`deposition_drive_Pa` 和 `grain_id`。
- `output/history.csv`：等效膜厚、最高膜面、partition-of-unity 误差和最大沉积驱动力。
- `output/run_parameters.csv`：本次运行的主要参数快照。

桥接层是必要的，因为当前 OpenPhase Python binding 暴露了 `DrivingForce` 的推进函数，却没有暴露 `Force(i,j,k).add_raw(...)`。相场求解本身仍由 OpenPhase C++ 库完成，Python 只负责算例编排和结果输出。

## 6. 使用方法

进入算例目录并运行快速检查：

```bash
cd /home/majie/trainFolder/openphaseFolder/02_ref_pvd/01_model_try
/home/majie/trainFolder/openphaseFolder/op_env/bin/python simulate.py \
  --steps 20 --output-every 10 --output smoke_output
```

第一次运行会自动执行 `make`。完整默认算例：

```bash
/home/majie/trainFolder/openphaseFolder/op_env/bin/python simulate.py
```

以相对法线 $35^\circ$ 的斜入射束流运行：

```bash
/home/majie/trainFolder/openphaseFolder/op_env/bin/python simulate.py \
  --incidence-angle 35 --flux 1.0 --flux-to-driving 2.0e7
```

完全关闭沉积可用于对照曲率松弛：

```bash
/home/majie/trainFolder/openphaseFolder/op_env/bin/python simulate.py \
  --flux 0 --steps 2000 --output no_deposition
```

用 ParaView 打开 `fields_*.vtk`，以 `grain_id` 着色观察晶粒竞争，以 `solid_fraction=0.5` 等值线观察膜面，以 `deposition_drive_Pa` 检查只有迎着通量的气固界面得到正驱动力。

## 7. 第一阶段应通过的检查

1. `history.csv` 中 `maximum_partition_error` 应接近浮点误差，不能随沉积系统增长。
2. 正常入射时等效膜厚应增长；`--flux 0` 对照组不应出现持续的净沉积增长。
3. 沉积驱动力只应出现在气固弥散界面，固相内部和晶界内部应为零。
4. 改变入射角后，`deposition_drive_Pa` 的高值区应转向迎风表面。
5. 把 `dt`、`dx` 或驱动力增大前，应分别进行时间步和网格收敛检查，并检查 OpenPhase 驱动力限制诊断。

## 8. 后续扩展顺序

第二阶段可增加 vapor density $\rho(\mathbf x,t)$ 的对流—扩散—消耗方程，并以 $\mathbf J=\rho\mathbf v$ 替代常通量；随后加入速度涨落和真正的 line-of-sight 遮挡。第三阶段再接入 SiMTRA 给出的角度、能量和通量分布，并加入材料特定的晶面各向异性、表面扩散、成核及实验参数标定。这样的分层实现可把 MPF/沉积耦合错误与输运或材料模型错误分开验证。
