# HelioSim: 优化的太阳能电池漂移扩散模型代码结构

## 核心文件

### 1. SolarCellParamsOptimized.m
- 太阳能电池参数类
- 包含物理常数、材料参数和网格设置
- 专门针对钙钛矿太阳能电池进行了优化
- 提供方便的参数设置方法

### 2. DDSolverChebfunOptimized.m
- 使用Chebfun谱方法的漂移扩散方程求解器
- 实现泊松方程和连续性方程的高精度求解
- 集成了界面处理和Scharfetter-Gummel离散化
- 使用自适应时间步长和Newton迭代

### 3. RecombinationModelsOptimized.m
- 优化的复合模型
- 实现SRH、Auger和辐射复合
- 针对界面复合进行了特殊处理
- 支持能量依赖的陷阱模型

### 4. OpticalGenerationOptimized.m
- 光生载流子生成模型
- 实现Beer-Lambert吸收
- 支持波长依赖的吸收系数
- 针对钙钛矿材料进行了优化

### 5. JVAnalyzerOptimized.m
- J-V特性曲线分析器
- 计算开路电压、短路电流、填充因子和效率
- 实现最大功率点跟踪
- 优化的电压扫描算法

### 6. VisualizerOptimized.m
- 结果可视化工具
- 能带图绘制
- 载流子密度可视化
- J-V曲线绘制
- 电场分布可视化

### 7. main_perovskite_cell.m
- 钙钛矿太阳能电池的主脚本
- 设置参数和配置
- 运行模拟
- 分析和可视化结果

## 数据流

1. `main_perovskite_cell.m` 创建 `SolarCellParamsOptimized` 实例并设置参数
2. 创建 `DDSolverChebfunOptimized` 实例，传入参数和配置
3. 求解器内部创建 `RecombinationModelsOptimized` 和 `OpticalGenerationOptimized` 实例
4. 求解器求解漂移扩散方程并返回结果
5. 创建 `JVAnalyzerOptimized` 实例分析性能
6. 创建 `VisualizerOptimized` 实例可视化结果

## 算法优化

1. **空间离散化**:
   - 使用Chebfun谱方法替代有限差分
   - 在界面处使用自适应网格细化
   - 实现高精度导数计算

2. **时间积分**:
   - 使用隐式时间步进
   - 实现自适应时间步长
   - 使用Newton方法求解非线性方程组

3. **界面处理**:
   - 使用Scharfetter-Gummel离散化
   - 考虑能带不连续性
   - 实现热电子发射模型

4. **复合模型**:
   - 实现位置依赖的复合参数
   - 界面特定复合
   - 陷阱能级分布

## 钙钛矿太阳能电池模型

### 结构
- ETL: TiO2 (100 nm)
- 吸收层: MAPbI3 (500 nm)
- HTL: Spiro-OMeTAD (100 nm)

### 关键参数
- **TiO2 (ETL)**:
  - 带隙: 3.2 eV
  - 电子亲和能: 4.0 eV
  - 介电常数: 9.0
  - 电子迁移率: 100 cm²/Vs
  - 施主掺杂: 1×10¹⁷ cm⁻³

- **MAPbI3 (吸收层)**:
  - 带隙: 1.55 eV
  - 电子亲和能: 3.9 eV
  - 介电常数: 25.0
  - 电子/空穴迁移率: 20 cm²/Vs
  - 本征(无掺杂)
  - 载流子寿命: 1 μs

- **Spiro-OMeTAD (HTL)**:
  - 带隙: 3.0 eV
  - 电子亲和能: 2.1 eV
  - 介电常数: 3.0
  - 空穴迁移率: 50 cm²/Vs
  - 受主掺杂: 1×10¹⁷ cm⁻³

### 界面参数
- ETL/吸收层界面复合速率: 10⁴ cm/s
- 吸收层/HTL界面复合速率: 10⁴ cm/s

## 预期性能
- 开路电压(Voc): ~1.0-1.1 V
- 短路电流(Jsc): ~22-24 mA/cm²
- 填充因子(FF): ~0.75-0.80
- 光电转换效率(PCE): ~18-22%
