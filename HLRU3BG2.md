---
title: "Flux-Balanced Patankar-type Schemes for the Compressible Euler Equations"
date: "2026-03-17T09:34:06+08:00"
draft: false
article_date: "2026-02-17"
authors:
  - "Thomas Izgin"
  - "Andreas Meister"
  - "Chi-Wang Shu"
  - "Davide Torlo"
categories:
  - "reading-notes"
tags:
  - "Mathematics - Numerical Analysis"
doi: "10.48550/arXiv.2602.14392"
url_source: "[http://arxiv.org/abs/2602.14392](http://arxiv.org/abs/2602.14392)"
zotero_key: "HLRU3BG2"
math: true
---
# 一、引言

流体模拟中，密度、内能等物理量必须保持正性，否则状态方程无意义。显式格式通常依赖 CFL 条件保证正性；要得到**无条件正性**（任意 $\Delta t$ 下保持正），格式必须是**非线性的**。Patankar 型方法通过引入非线性权重实现这一点，同时保持守恒性。MPRK（Modified Patankar–Runge–Kutta）即基于 Runge–Kutta 框架的 Patankar 型格式。

---

# 二、PDRS 理论基础

## 2.1 定义

MPRK 针对 **Production-Destruction-Rest System (PDRS)** 设计，形式为

$$
u'_i(t) = r_i(u,t) + \sum_{j=1}^{N} \bigl( p_{ij}(u,t) - d_{ij}(u,t) \bigr), \quad u(0) = u_0 \in \mathbb{R}^N_{>0}

$$

其中：

- $p_{ij}, d_{ij} \geq 0$：产生项与消灭项
- $r_i = r^p_i - r^d_i$，$r^p_i, r^d_i \geq 0$：剩余项
- **守恒配对**：$p_{ij} = d_{ji}$，$p_{ii} = d_{ii} = 0$

## 2.2 通量分解：从半离散到 PDRS

对半离散形式 $\frac{du_i}{dt} = \frac{1}{\Delta x}(\hat{F}_{i-1/2} - \hat{F}_{i+1/2})$，利用恒等式

$$
\hat{F} = \max\{0,\hat{F}\} - \bigl(-\min\{0,\hat{F}\}\bigr)

$$

将通量差写成非负产生减非负消灭。界面 $i+\frac{1}{2}$ 上：

- $\hat{F}_{i+1/2} > 0$（左→右）：$d_{i,i+1} = p_{i+1,i} = \hat{F}_{i+1/2}$
- $\hat{F}_{i+1/2} < 0$（右→左）：$p_{i,i+1} = d_{i+1,i} = -\hat{F}_{i+1/2}$

边界通量同理拆成 $r^p$、$r^d$。流出单元 $i$ 的贡献等于流入相邻单元的贡献，满足 $p_{ij} = d_{ji}$。

---

# 三、MPRK 格式

## 3.1 格式定义

给定 s 级显式 RK 方法（Butcher 表 **$A,b,c \geq 0$**），对应的 MPRK 格式为

**中间步**（$k = 1,\ldots,s$）：

$$
u^{(k)}_i = u^n_i + \Delta t \sum_{\nu=1}^{k-1} a_{k\nu} \Biggl(
  r^p_i + \sum_j p_{ij} \frac{u^{(k)}_j}{\pi^{(k)}_j}

- \Bigl( r^d_i + \sum_j d_{ij} \Bigr) \frac{u^{(k)}_i}{\pi^{(k)}_i}
\Biggr)

$$

**最终步**：

$$
u^{n+1}_i = u^n_i + \Delta t \sum_{k=1}^{s} b_k \Biggl(
  r^p_i + \sum_j p_{ij} \frac{u^{n+1}_j}{\sigma_j}

- \Bigl( r^d_i + \sum_j d_{ij} \Bigr) \frac{u^{n+1}_i}{\sigma_i}
\Biggr)

$$

$\pi^{(k)}_i$、$\sigma_i$ 为 **Patankar 权分母 (PWDs)**，须 $\pi^{(k)}_i, \sigma_i > 0$。产生项乘 $u^{n+1}_j/\sigma_j$，消灭项乘 $u^{n+1}_i/\sigma_i$。

## 3.2 PWD 选取


| 格式               | PWD 选取                                      | 说明           |
| -------------------- | ----------------------------------------------- | ---------------- |
| **MPE**（一阶）    | $\sigma_i = u^n_i$                            | 基于前向 Euler |
| **MPHeun**（二阶） | $\pi^{(2)}_i = u^n_i$，$\sigma_i = u^{(2)}_i$ | 基于 Heun 法   |

Patankar 权 $\frac{u^{n+1}_i}{u^n_i}$ 在 $\Delta t \to 0$ 时趋于 1，故不破坏收敛阶。

## 3.3 矩阵形式与 $M^{(k)}$ 的构成

移项后可写成 $M^{(k)} u^{(k)} = u^n + \Delta t \sum_\nu a_{k\nu} r^p$ 及 $M u^{n+1} = u^n + \Delta t \sum_k b_k r^p$。

**$M^{(k)}$ 的构成**：将右端所有乘在 $u^{(k)}$ 上的项移到左端，左端 $u^{(k)}$ 的系数即 $M^{(k)}$ 的元素。对角元来自 $1$ 与消灭项，非对角元来自产生项（带负号）。

## 3.4 格式性质


| 性质             | 说明                                                                                                                                 |
| ------------------ | -------------------------------------------------------------------------------------------------------------------------------------- |
| **隐格式**       | 未知量$u^{(k)}$、$u^{n+1}$ 在等式两端，需解线性系统；为**线性隐式**（无需 Newton 迭代）                                              |
| **M-矩阵**       | 当$M^{(k)}$、$M$ 为 M-矩阵时无条件正性。$\Delta t$ 充分小时 $M \approx I$ 自然成立；任意 $\Delta t$ 下需 PWD 与 $p_{ij}=d_{ji}$ 保证 |
| **Butcher 非负** | 要求$a_{k\nu}, b_k, c_k \geq 0$，否则破坏 $m_{ij} \leq 0$，无法保证 M-矩阵。适用 Euler、Heun、SSP-RK 等                              |

## 3.5 从 PDRS 到 MPE

取 $\sigma_i = u^n_i$，MPE 的更新为

$$
u^{n+1}_i = u^n_i + \Delta t r^p_i

- \Delta t \left( \sum_j p_{ij} \frac{u^{n+1}_j}{u^n_j}

- \Bigl(r^d_i + \sum_j d_{ij}\Bigr) \frac{u^{n+1}_i}{u^n_i} \right)

$$

移项得 $M u^{n+1} = u^n + \Delta t r^p$，$M$ 为三对角 M-矩阵。

---

# 四、密度方程的 MPE 应用

## 4.1 半离散方程 (16)

$$
\frac{d\rho_i}{dt} = \frac{1}{\Delta x}\left(\hat{F}_{i-1/2,1} - \hat{F}_{i+1/2,1}\right), \quad i = 1, \ldots, N

$$

$S_1(U_i) = 0$。非周期边界下 $\sum_i (\hat{F}_{i-1/2,1} - \hat{F}_{i+1/2,1}) \neq 0$，故需写成 PDRS。

## 4.2 PDRS 形式 (17)

$$
\frac{d\rho_i}{dt} = \frac{1}{\Delta x}\left(r^p_{i,1} - r^d_{i,1} + \sum_{j=1}^{N}(p_{i,j,1} - d_{i,j,1})\right)

$$

其中 $r^p_{i,1}, r^d_{i,1}, p_{i,j,1}, d_{i,j,1} \geq 0$，$p_{i,j,1} = d_{j,i,1}$：

- 内界面：$p_{i+1,i,1} = d_{i,i+1,1} = \max\{0, \hat{F}_{i+1/2,1}\}$，$p_{i,i+1,1} = d_{i+1,i,1} = -\min\{0, \hat{F}_{i+1/2,1}\}$
- $|i-j| \neq 1$ 时 $p_{i,j,1} = d_{j,i,1} = 0$

边界项 (19)：

$$
r^p_1 = \max\{0, \hat{F}_{1/2,1}\}\, e_1 + \max\{0, -\hat{F}_{N+1/2,1}\}\, e_N

$$

$$
r^d_1 = -\left(\min\{0, \hat{F}_{1/2,1}\}\, e_1 + \min\{0, -\hat{F}_{N+1/2,1}\}\, e_N\right)

$$

## 4.3 $i=1$ 时的具体形式

$$
\begin{aligned}
r^p_{1,1} &= \max\{0,\, \hat{F}_{\frac{1}{2},1}\}, \quad
r^d_{1,1} = -\min\{0,\, \hat{F}_{\frac{1}{2},1}\} \\[4pt]
p_{1,j,1} &= \begin{cases} -\min\{0,\, \hat{F}_{\frac{3}{2},1}\}, & j=2 \\ 0, & j\neq 2 \end{cases} \\[4pt]
d_{1,j,1} &= \begin{cases} \max\{0,\, \hat{F}_{\frac{3}{2},1}\}, & j=2 \\ 0, & j\neq 2 \end{cases}
\end{aligned}

$$

## 4.4 MPE 更新 (20)

取 $\sigma_i = \rho^n_i$，

$$
\rho^{n+1}_i = \rho^n_i + \frac{\Delta t}{\Delta x} r^{p,n}_{i,1}
+ \frac{\Delta t}{\Delta x}\left(\sum_{j=1}^{N} p^n_{i,j,1}\,\frac{\rho^{n+1}_j}{\rho^n_j}
- \Bigl(r^{d,n}_{i,1} + \sum_{j=1}^{N} d^n_{i,j,1}\Bigr)\frac{\rho^{n+1}_i}{\rho^n_i}\right)

$$

产生项乘 $\rho^{n+1}_j/\rho^n_j$，消灭项乘 $\rho^{n+1}_i/\rho^n_i$，保证 unconditionally positive。

---

# 五、源项的 Patankar 加权（反应 Euler）

对耦合反应 $S_k = S^p_k - S^d_k$，在 MPE 中加源项时采用 Patankar 权：

**$U_{i,1}$ 方程**：

- $S^p_1$ 乘 $\Delta t\frac{U_{i,2}^{n+1}}{U_{i,2}^n}$（产生来自 $U_2$）
- $S^d_1$ 乘 $\Delta t\frac{U_{i,1}^{n+1}}{U_{i,1}^n}$（消灭针对 $U_1$）

**$U_{i,2}$ 方程**：

- $S^p_2$ 乘 $\Delta t\frac{U_{i,1}^{n+1}}{U_{i,1}^n}$（产生来自 $U_1$）
- $S^d_2$ 乘 $\Delta t\frac{U_{i,2}^{n+1}}{U_{i,2}^n}$（消灭针对 $U_2$）

守恒关系 $S^p_1 = S^d_2$、$S^d_1 = S^p_2$ 使两方程源项之和为零，保持质量守恒。权重 $\frac{U^{n+1}}{U^n}$ 实现交叉耦合与自反馈，有助于 stiff 源项的数值稳定。

---

# 六、MP trick、盲目应用之弊端与 flux-balanced 策略（小结）

## MP trick 与无条件正性

**MP trick（Patankar trick）**：在产生项乘 $u_j^{n+1}/u_j^n$、消灭项乘 $u_i^{n+1}/u_i^n$，使格式变为隐式，且 $M$ 为 M-矩阵，从而在任意 $\Delta t > 0$ 下保证 $u^{n+1} > 0$。代价是格式非线性。

## 盲目应用的 drawback

若对**每个正守恒变量**（密度 $\rho$ 与总能量 $\rho E$）都使用 MP trick（即 MPE-ρE、MPHeun-ρE）：

- **理论**：为保持 contact discontinuity（$u$、$p$ 常数），动量、能量通量须与密度通量采用**同一组** Patankar 权重；而对能量单独做 MP trick 会引入**另一套**权重，二者冲突。
- **数值**：破坏接触间断（Figure 2）；真空算例出现非物理振荡与负压（Table 5）；CFL 限制更严、效率更低（Table 3、4）。

## Flux-balanced 策略：仅对密度做 MP trick

文中的最终方案：**只对密度方程用 MP trick**，动量、能量通量**不**单独做 MP trick，而是按**密度方程的 Patankar 权重** $w_i = \rho^{n+1}_i / \rho^n_i$ 对相关通量部分加权。

- 能量通量拆成 $\hat{F}^w + \hat{F}^u$：与密度通量结构类似的部分用 $w_i$ 加权，其余保持显式。
- **$\rho E$ 的处理**：不做 MP trick，无理论上的无条件正性；但数值上压强正性更好、接触间断保持、CFL 更大（如 MPHeun-s 可达 CFL ≈ 0.86），且 momentum/energy 不再额外解线性系统。

---

# 七、CR格式与PDRS的结合

界面、边界和源项也许可以考虑PDRS，产生-消耗-残留。

## 6.1 分解结构

CR 高阶格式的 antidiffusive flux 可写为

$$
f_i = f^K_i + f^F_i, \quad
f^K_i = \sum_{K} \int_K (f(u_h) - f_h) \cdot \nabla\varphi_i \, dx, \quad
f^F_i = -\sum_{F} \int_F \hat{f}(u_h^-, u_h^+; n_F) \llbracket\varphi_i\rrbracket \, ds

$$

单元部分 $f^K_i$ 保持原有 FCT：$f^K_{ij} = \frac{1}{3}(f^K_i - f^K_j)$，$f^K_{ij} = -f^K_{ji}$，再对 $g^K_{ij} = f^K_{ij} - d^n_{ij}(U^n_j - U^n_i)$ 做 limiter。下面**仅对面通量部分** $f^F_i$ 设计 PDRS 算法。

## 6.2 面通量的 PDRS 分解

对每条内边 $F \in \mathcal{F}_h^i$，记 $\hat{f}_F := \int_F \hat{f}(u_h^-, u_h^+; n_F)\, ds$（或中点求积）。$n_F$ 由 $F^-$ 指向 $F^+$。记 $I^-(F)$ 为 $F^-$ 侧单元的边指标集，$I^+(F)$ 为 $F^+$ 侧边指标集，满足 $I(F) = I^-(F) \cup I^+(F)$，且公共边 $F$ 可任归一侧。

**步骤 1–2**：将 $\hat{f}_F$ 分配到 DOF 对。设权重 $w^+_i,\, w^-_j \geq 0$，$\sum_{i \in I^+(F)} w^+_i = \sum_{j \in I^-(F)} w^-_j = 1$。定义配对通量（$j \in I^-(F) \to i \in I^+(F)$ 表示从 $F^-$ 到 $F^+$）

$$
f^F_{ij} := \hat{f}_F \cdot w^+_i \cdot w^-_j \quad (i \in I^+(F),\, j \in I^-(F)), \qquad
f^F_{ji} := -f^F_{ij}

$$

再按 Izgin 做符号分解：$p^F_{ij} := \max\{0, f^F_{ij}\}$，$d^F_{ij} := \max\{0, -f^F_{ij}\}$。则 $p^F_{ij} - d^F_{ij} = f^F_{ij}$ 且 $p^F_{ij} = d^F_{ji}$。

**权重选取**：为与 $f^F_i = -\int_F \hat{f}\llbracket\varphi_i\rrbracket$ 一致，可取

$$
w^{\pm}_i := \frac{m_i \, |\llbracket\varphi_i\rrbracket_F|}{\sum_{k \in I^{\pm}(F)} m_k \, |\llbracket\varphi_k\rrbracket_F|}

$$

其中 $\llbracket\varphi_i\rrbracket_F$ 为 $\varphi_i$ 在 $F$ 上的跳变（一侧为 0 时取另一侧绝对值）。分母为 0 时改用均匀权重。

## 6.3 时间步格式与完整算法

对 DOF $i$，记 $J^+_i := \{ (F,j) : F \ni i,\, i \in I^+(F),\, j \in I^-(F) \}$，$J^-_i := \{ (F,j) : F \ni i,\, i \in I^-(F),\, j \in I^+(F) \}$。面通量的 Patankar 修正为

$$
C^F_i = \frac{\Delta t}{m_i} \Biggl(
  \sum_{(F,j) \in J^+_i} p^F_{ij} \frac{U^{n+1}_j}{U^n_j}
  + \sum_{(F,j) \in J^-_i} p^F_{ji} \frac{U^{n+1}_i}{U^n_i}
  - \sum_{(F,j) \in J^-_i} d^F_{ij} \frac{U^{n+1}_i}{U^n_i}
  - \sum_{(F,j) \in J^+_i} d^F_{ji} \frac{U^{n+1}_j}{U^n_j}
\Biggr)

$$

由 $p^F_{ij} = d^F_{ji}$，可统一写成（与 Izgin 一致）

$$
C^F_i = \frac{\Delta t}{m_i} \Biggl( \sum_j p^F_{ij} \frac{U^{n+1}_j}{U^n_j}
- \Bigl( r^{d,F}_i + \sum_j d^F_{ij} \Bigr) \frac{U^{n+1}_i}{U^n_i} \Biggr)

$$

其中 $r^{d,F}_i$ 为边界面贡献，$p^F_{ij},\, d^F_{ij}$ 由步骤 1–2 组装。完整更新为

$$
U^{n+1}_i = U^{n+1,L}_i + \frac{\Delta t}{m_i} \sum_j l^n_{ij} g^n_{ij} + C^F_i

$$

**算法流程**：

1. 计算低阶解 $U^{n+1,L}$；
2. 计算 $f^K_i$，单元分解 $f^K_{ij}$，形成 $g^K_{ij} = f^K_{ij} - d^n_{ij}(U^n_j - U^n_i)$，应用 limiter 得 $l^K_{ij}$；
3. 对每条内边 $F$：计算 $\hat{f}_F$，按权重得 $f^F_{ij}$，再得 $p^F_{ij} = \max\{0, f^F_{ij}\}$、$d^F_{ij} = \max\{0, -f^F_{ij}\}$；
4. 求解耦合系统：$U^{n+1}_i = U^{n+1,L}_i + \frac{\Delta t}{m_i}\sum_j l^K_{ij} g^K_{ij} + C^F_i$，其中 $C^F_i$ 含 $U^{n+1}/U^n$，为线性隐式（每次迭代或直接求解 $M U^{n+1} = b$）。

### 6.3.1 $l_{ij}$ 的保正构造

将完整更新移项后写成 $M U^{n+1} = \tilde{b}$，其中矩阵 $M$ 来自 $C^F_i$（定理 6.1 已证为 M-矩阵），右端为

$$
\tilde{b}_i = U^{n+1,L}_i + \frac{\Delta t}{m_i} \sum_j l_{ij} g_{ij}.
$$

因 $M^{-1} \geq 0$，要 $U^{n+1} = M^{-1} \tilde{b} \geq 0$，只需 **$\tilde{b} \geq 0$**。故 $l_{ij}$ 须满足

$$
U^{n+1,L}_i + \frac{\Delta t}{m_i} \sum_j l_{ij} g_{ij} \geq 0 \quad \forall i.
$$

由 $g_{ij} = -g_{ji}$ 知单元修正守恒。对节点 $i$，当 $g_{ij} < 0$ 时，$(Δt/m_i) l_{ij} g_{ij} < 0$ 会使 $\tilde{b}_i$ 下降。故**关键约束**为：流出 $i$ 的反扩散通量之和不能超过 $U^{n+1,L}_i$：

$$
\frac{\Delta t}{m_i} \sum_{j \,:\, g_{ij} < 0} l_{ij} |g_{ij}| \leq U^{n+1,L}_i.
$$

**构造思路**（对称 $l_{ij} = l_{ji} \in [0,1]$）：下面三种方式均保证 $\tilde{b} \geq 0$ 且守恒。

---

**方式一：单点比例法**

对每个节点 $i$，记其“流出容量”（$l=1$ 时的最大可能流出）
$$
Q^-_i := \sum_{j \,:\, g_{ij} < 0} |g_{ij}|.
$$
保正约束 $(Δt/m_i) \sum_{j:g_{ij}<0} l_{ij}|g_{ij}| \leq U^{n+1,L}_i$ 等价于：若对所有 $j \sim i$ 使用同一因子 $\theta_i$，则需 $\theta_i Q^-_i \leq U^{n+1,L}_i m_i/\Delta t$。定义
$$
\theta_i := \begin{cases}
\min\bigl\{1,\, \dfrac{U^{n+1,L}_i \, m_i}{\Delta t \, Q^-_i}\bigr\}, & Q^-_i > 0 \\
1, & Q^-_i = 0
\end{cases}
$$
对每条边 $(i,j)$：当 $g_{ij} < 0$ 时，$i$ 为捐出端，限制来自 $\theta_i$；当 $g_{ij} > 0$ 时，$j$ 为捐出端，限制来自 $\theta_j$。取对称
$$
l_{ij} = l_{ji} := \min\{\theta_i,\, \theta_j\}
$$
即同时满足 $i,j$ 的约束，且 $l_{ij} = l_{ji}$ 保持守恒。验证：对 $i$，$\sum_{j:g_{ij}<0} l_{ij}|g_{ij}| \leq \theta_i Q^-_i \leq U^{n+1,L}_i m_i/\Delta t$。

---

**方式二：逐对 Zalesak 型（推荐）**

对每条边 $(i,j)$ 单独限幅。记节点 $k$ 可捐出量为 $R^-_k := U^{n+1,L}_k$。边 $(i,j)$ 上的反扩散通量为 $(Δt/m_i) l_{ij} g_{ij}$（对 $i$）与 $(Δt/m_j) l_{ji} g_{ji} = -(Δt/m_j) l_{ij} g_{ij}$（对 $j$，因 $g_{ji}=-g_{ij}$）。

- **$g_{ij} > 0$**：通量 $j \to i$。$j$ 捐出 $(Δt/m_j) l_{ij} g_{ij}$，$i$ 接收。约束：$(Δt/m_j) l_{ij} g_{ij} \leq R^-_j$，故
  $$
  l_{ij} \leq \dfrac{R^-_j \, m_j}{\Delta t \, g_{ij}} = \dfrac{U^{n+1,L}_j \, m_j}{\Delta t \, g_{ij}}.
  $$
- **$g_{ij} < 0$**：通量 $i \to j$。$i$ 捐出 $(Δt/m_i) l_{ij} |g_{ij}|$。约束：$(Δt/m_i) l_{ij} |g_{ij}| \leq R^-_i$，故
  $$
  l_{ij} \leq \dfrac{R^-_i \, m_i}{\Delta t \, |g_{ij}|} = \dfrac{U^{n+1,L}_i \, m_i}{\Delta t \, |g_{ij}|}.
  $$

综合并保持 $l_{ij} \in [0,1]$：
$$
l_{ij} = l_{ji} := \begin{cases}
\min\Bigl\{1,\, \dfrac{U^{n+1,L}_j \, m_j}{\Delta t \, g_{ij}}\Bigr\}, & g_{ij} > 0 \\[6pt]
\min\Bigl\{1,\, \dfrac{U^{n+1,L}_i \, m_i}{\Delta t \, |g_{ij}|}\Bigr\}, & g_{ij} < 0 \\
1, & g_{ij} = 0
\end{cases}
$$
（$g_{ij}=0$ 时该边无贡献，取 $l_{ij}=1$ 无影响。）

---

**方式三：迭代裁剪**

初值：$l_{ij} = l_{ji} = 1$ 对所有边 $(i,j)$。

1. 计算 $\tilde{b}_i = U^{n+1,L}_i + \frac{\Delta t}{m_i} \sum_j l_{ij} g_{ij}$；
2. 若 $\tilde{b}_k \geq 0$ 对所有 $k$，停止；
3. 否则取 $k$ 使得 $\tilde{b}_k < 0$。记 $k$ 的流入、流出为
   $$
   P_k := \sum_{j \,:\, g_{kj} > 0} l_{kj} g_{kj}, \qquad
   D_k := \sum_{j \,:\, g_{kj} < 0} l_{kj} |g_{kj}|.
   $$
   则 $\tilde{b}_k = U^{n+1,L}_k + \frac{\Delta t}{m_k}(P_k - D_k)$。若将流出端系数 $l_{kj}$（$g_{kj}<0$）统一乘以 $\theta_k$，则新流出为 $\theta_k D_k$，新右端为
   $$
   \tilde{b}_k^{\mathrm{new}} = U^{n+1,L}_k + \frac{\Delta t}{m_k}(P_k - \theta_k D_k).
   $$
   令其 $\geq 0$ 得 $\theta_k D_k \leq U^{n+1,L}_k m_k/\Delta t + P_k$，即
   $$
   \theta_k := \min\Bigl\{1,\, \frac{U^{n+1,L}_k \, m_k/\Delta t + P_k}{D_k}\Bigr\}
   $$
   （$D_k=0$ 时取 $\theta_k=1$）。对所有 $j \in \mathcal{N}^-_k := \{j : g_{kj} < 0\}$ 做 $l_{kj} \leftarrow \theta_k \, l_{kj}$，$l_{jk} \leftarrow \theta_k \, l_{jk}$。返回步骤 1。

**说明**：缩减 $l_{kj}$ 会减少 $k$ 的流出（提升 $\tilde{b}_k$），但同时也减少 $j$ 的流入（可能使 $\tilde{b}_j$ 降低）。故单次修正 $k$ 后，可能与 $k$ 相邻的 $j$ 变为负；需多轮遍历直至全体 $\tilde{b} \geq 0$。实际中通常数轮即可。

---

> ex_article 的 monolithic convex limiting（Remark 4.2）在无 PDRS 时保 DMP；此处需额外保证 $\tilde{b} \geq 0$。在 $U^n > 0$ 且 CFL 满足时 $U^{n+1,L} > 0$，上述三种方式均可达成 $\tilde{b} \geq 0$。

### 6.3.2 向保极值的扩展

记 $U_{\min} := \min_j U^n_j$，$U_{\max} := \max_j U^n_j$。要得到 **DMP**（$U^{n+1}_i \in [U_{\min},\, U_{\max}]$），有三条可行路径。

---

**思路 A：加强 $l_{ij}$ 约束，使 $\tilde{b} \in [U_{\min},\, U_{\max}]$**

将保正条件 $U_{\min} \leq \tilde{b}_i$（当 $U_{\min} > 0$ 时比 $\tilde{b}_i \geq 0$ 更紧）和上界 $\tilde{b}_i \leq U_{\max}$ 一并施加：

$$
U_{\min} \leq U^{n+1,L}_i + \frac{\Delta t}{m_i} \sum_j l_{ij} g_{ij} \leq U_{\max} \quad \forall i.
$$

记 $P_i = \sum_{j:g_{ij}>0} l_{ij} g_{ij}$，$D_i = \sum_{j:g_{ij}<0} l_{ij}|g_{ij}|$，则：

- **下界**：$(Δt/m_i) D_i \leq U^{n+1,L}_i - U_{\min} + (Δt/m_i) P_i$
- **上界**：$(Δt/m_i)(P_i - D_i) \leq U_{\max} - U^{n+1,L}_i$

以**方式二**为例，保极值扩展为：边 $(i,j)$ 上捐出端最多捐出到边界的距离，接收端最多接收到边界的距离。取

- $g_{ij} > 0$（$j \to i$）：捐出端 $j$ 可捐 $\leq U^{n+1,L}_j - U_{\min}$，接收端 $i$ 可收 $\leq U_{\max} - U^{n+1,L}_i$，故
  $$
  l_{ij} = \min\Bigl\{1,\, \frac{(U^{n+1,L}_j - U_{\min}) m_j}{\Delta t \, g_{ij}},\, \frac{(U_{\max} - U^{n+1,L}_i) m_i}{\Delta t \, g_{ij}}\Bigr\};
  $$
- $g_{ij} < 0$（$i \to j$）：对称地
  $$
  l_{ij} = \min\Bigl\{1,\, \frac{(U^{n+1,L}_i - U_{\min}) m_i}{\Delta t \, |g_{ij}|},\, \frac{(U_{\max} - U^{n+1,L}_j) m_j}{\Delta t \, |g_{ij}|}\Bigr\}.
  $$

方式一、三可类似将 $R^-_k$ 从 $U^{n+1,L}_k$ 改为 $\min\{U^{n+1,L}_k - U_{\min},\, U_{\max} - U^{n+1,L}_k\}$，并加上接收端约束。

**注意**：$\tilde{b} \in [U_{\min},\, U_{\max}]$ 只能控制右端；最终解为 $U^{n+1} = M^{-1} \tilde{b}$。当 $M$ 为 M-矩阵且**行和为 1**（$M \mathbf{1} = \mathbf{1}$）时，$M^{-1}$ 非负且每行和为 1，故凸组合保持区间，从而 $U^{n+1} \in [U_{\min},\, U_{\max}]$。PDRS 矩阵 $M$ 是否满足行和 1，需根据 $p^F_{ij} = d^F_{ji}$ 及守恒性具体验证；若不满足，则需思路 B。

---

**思路 B：后处理裁剪（6.4 已给）**

$$
U^{n+1}_i \leftarrow \operatorname{median}\bigl(U^{n+1}_i,\, U_{\min},\, U_{\max}\bigr)
$$

实现简单，但会破坏守恒（裁剪会改变 $\sum_i m_i U^{n+1}_i$）。若守恒重要，可改为：先裁剪，再对超限节点做局部守恒修正（如 scaling 或 flux redistribution），使总质量不变。

---

**思路 C：$C^F$ 较小时的近似**

当面通量 $C^F_i$ 相对单元修正较小时，$M \approx I$，$U^{n+1} \approx \tilde{b}$。此时只要 $\tilde{b} \in [U_{\min},\, U_{\max}]$，即有 $U^{n+1} \in [U_{\min},\, U_{\max}]$，思路 A 的 $l_{ij}$ 约束在实践上足够。

---

**小结**：若 PDRS 矩阵行和为 1，可将三种方式的 $R^-_k$ 改为到边界的距离，得到保极值格式；否则建议用后处理 median，并视需要做守恒修正。

