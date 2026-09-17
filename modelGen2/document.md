# OpenPhase 多晶 PVD + vapor density 对流—扩散模型

1. 模型范围
-----------
该程序在第一阶段 vapor->grain 成对沉积驱动力的基础上，增加独立的归一化
vapor density rho(x,t)。OpenPhase 仍负责 Steinbach 双障碍势多相场、界面能、
晶界和相分数约束；Python 使用显式有限体积格式推进 vapor transport。

当前模型包含：多晶粒相场、定向沉积、vapor 扩散、可选对流、沉积消耗和顶部
持续供料。尚未包含 SiMTRA、真实原子数密度标定、表面扩散、成核、弹性和射线
追踪意义下的全局遮挡。

2. 控制方程
-----------
令 phi_v 为 OpenPhase 聚合气相分数，phi_s=sum_i(phi_i) 为总固相分数，满足

    phi_v + sum_i(phi_i) = 1.

脚本求解正则化的掩膜对流—扩散—消耗方程

    d rho/dt = div(D phi_v grad(rho))
               - div(phi_v rho v)
               - C_dep max(d phi_s/dt, 0).

x 方向周期；顶部 rho=rho_feed；底部无通量。面通量乘相邻单元 phi_v 的较小值，
防止 vapor transport 穿过致密固相。对流关闭时令 v=0，扩散和沉积仍然计算。

局部薄膜外法向为

    n = grad(phi_v)/|grad(phi_v)|,

入射方向因子为

    q = max(-Jhat dot n, 0).

写入 OpenPhase 的 vapor->grain 成对沉积驱动力为

    DeltaG_vi = Lambda * rho * |J| * q * w_i,
    w_i = phi_i^p / sum_j(phi_j^p).

反对称 pairwise 更新使气相减少量等于固相增加量，不破坏相分数和为 1 的约束。

3. 数值方法和稳定性
-------------------
扩散项使用中心差分，对流项使用一阶迎风面通量。脚本根据

    4 D dt/dx^2 + (|vx|+|vy|) dt/dx

自动选择 transport 子步数，使每个子步的指标不超过 0.8。默认 dx=2 nm、
dt=1e-8 s、D=2e-11 m^2/s、|v|=0.02 m/s，通常只需一个 transport 子步。
一阶迎风具有数值扩散；后续若研究尖锐 vapor front，可换成 TVD/WENO 格式。

rho 当前是归一化密度，rho_feed 默认等于 1。consumption-factor 把固相分数增量
映射为 rho 消耗，它与 flux-to-driving 都需要用实际沉积速率或原子体积标定，
因此当前输出用于方法验证和参数趋势研究，不应直接解释成定量 PVD 工艺预测。

4. 文件
-------
simulate.py          主程序及 vapor transport 求解器。
ProjectInput.opi     OpenPhase 网格、相、界面能、迁移率和边界条件。
openphase_vapor.cpp  将 rho 依赖的沉积驱动力写入 OpenPhase DrivingForce。
Makefile             编译 pybind11 桥接模块。
document.py          本说明；执行后在终端打印文档。

输出 fields_*.vtk 包含 vapor_fraction、solid_fraction、vapor_density、
deposition_drive_Pa 和 grain_id。history.csv 记录膜厚、平均/总 vapor density、
累计沉积消耗、相分数约束误差及 transport 子步数。run_parameters.csv 保存完整
运行参数。

5. 使用方法
-----------
进入目录：

    cd /home/majie/trainFolder/openphaseFolder/02_ref_pvd/02_vapor_density

快速检查（默认开启对流）：

    /home/majie/trainFolder/openphaseFolder/op_env/bin/python simulate.py \
        --steps 20 --output-every 10 --output smoke_convection

关闭对流、只保留扩散：

    /home/majie/trainFolder/openphaseFolder/op_env/bin/python simulate.py \
        --no-convection --steps 20 --output-every 10 \
        --output smoke_diffusion_only

斜入射并改变 transport 参数：

    /home/majie/trainFolder/openphaseFolder/op_env/bin/python simulate.py \
        --incidence-angle 35 --transport-speed 0.03 \
        --vapor-diffusivity 1.0e-11

查看本说明：

    less document.md

6. 建议检查
-----------
有/无对流算例应均保持 maximum_partition_error 接近机器精度；vapor_density 不得
出现负值或超过 rho_feed；沉积驱动力只应出现在迎流气固界面；关闭 flux 后不应
出现持续净沉积。定量研究前应分别进行 dx、dt、transport 子步、D、速度和驱动力
的收敛性与敏感性检查。
