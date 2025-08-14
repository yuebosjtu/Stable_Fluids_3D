# 3D Stable Fluids Simulation

这是一个基于 Jos Stam 的 Stable Fluids 方法的三维流体模拟实现。

## 特性

- **三维 MAC 网格**: 使用交错网格(MAC grid)存储速度分量
- **半拉格朗日平流**: 实现无条件稳定的平流算法
- **压力投影**: 使用 AMG 预条件共轭梯度法求解压力
- **球形源**: 支持在三维空间中的球形速度源
- **边界条件**: 实现无滑移壁面边界条件
- **VTK 输出**: 生成可在 ParaView 中可视化的 VTK 文件

## 编译和运行

### 环境要求

- CMake 3.16 或更高版本
- C++17 兼容的编译器
- Eigen3 库（已包含在项目中）

### 编译步骤

```bash
mkdir build
cd build
cmake ..
cmake --build . --config Release
```

### 运行模拟

```bash
./main
```

## 项目结构

```
include/
├── stable_fluid.h      # 三维流体算法核心函数
├── fluid_simulator.h   # 主模拟器类定义
├── amg_solver.h       # AMG求解器（压力求解）
├── utils.h            # 工具函数（源场创建等）
└── Eigen/             # Eigen库头文件

src/
└── fluid_simulator.cpp # 主模拟器类实现

main.cpp               # 主程序入口
CMakeList.txt          # CMake构建文件
```

## 算法说明

### 1. MAC 网格布局

- **u 速度分量**: 存储在 `(i, j+0.5, k+0.5)` 位置，大小 `(width+1) × height × depth`
- **v 速度分量**: 存储在 `(i+0.5, j, k+0.5)` 位置，大小 `width × (height+1) × depth`  
- **w 速度分量**: 存储在 `(i+0.5, j+0.5, k)` 位置，大小 `width × height × (depth+1)`
- **压力和染料**: 存储在单元中心 `(i+0.5, j+0.5, k+0.5)`，大小 `width × height × depth`

### 2. 模拟步骤

每个时间步包含以下步骤：

1. **施加外力**: 重力等外部力
2. **速度平流**: 使用半拉格朗日方法平流速度场
3. **压力求解**: 求解泊松方程使速度场无散
4. **染料平流**: 平流染料/示踪场
5. **边界条件**: 应用壁面无滑移条件
6. **源条件**: 在球形区域内维持恒定速度

### 3. 三线性插值

实现了针对 MAC 网格的三线性插值：
- `interpolate_u()`: u分量插值
- `interpolate_v()`: v分量插值  
- `interpolate_w()`: w分量插值
- `trilerp_scalar()`: 标量场插值

### 4. 压力求解

使用 AMG 预条件共轭梯度法求解离散化的泊松方程：
```
∇²p = -∇·u
```

## 可视化

模拟生成以下输出：

1. **文本数据文件**:
   - `velocity_data.txt`: 速度场数据
   - `pressure_data.txt`: 压力场数据  
   - `dye_data.txt`: 染料浓度数据

2. **VTK 文件**: 
   - `velocity_*.vtk`: 可在 ParaView 中打开的三维可视化文件

## 参数调整

在 `main.cpp` 中可以调整以下参数：

- `width, height, depth`: 网格分辨率
- `dt`: 时间步长
- `total_time`: 总模拟时间
- `source_center_*`: 球形源位置
- `source_radius`: 球形源半径
- `source_velocity_*`: 源区域速度
- `pressure_iterations`: 压力求解迭代次数

## 性能建议

- 对于初次测试，建议使用较低分辨率（如 32³）
- 高分辨率模拟（如 128³）需要大量内存和计算时间
- 可以通过减少压力求解迭代次数来加速（可能影响精度）

## 扩展功能

未来可以添加的功能：
- MacCormack 平流方法
- 粘性扩散
- 温度对流
- 多相流体
- GPU 加速

## 许可证

本项目基于原始的 Stable Fluids 算法实现，仅供学习和研究使用。
