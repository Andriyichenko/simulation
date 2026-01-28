# Simulation Project Overview

![VS Code](https://img.shields.io/badge/Visual%20Studio%20Code-0078d7.svg?style=for-the-badge&logo=visual-studio-code&logoColor=white) ![C++](https://img.shields.io/badge/c++-%2300599C.svg?style=for-the-badge&logo=c%2B%2B&logoColor=white) ![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54) ![Jupyter Notebook](https://img.shields.io/badge/jupyter-%23FA0F00.svg?style=for-the-badge&logo=jupyter&logoColor=white) ![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white) ![Markdown](https://img.shields.io/badge/markdown-%23000000.svg?style=for-the-badge&logo=markdown&logoColor=white) ![Code Runner](https://img.shields.io/badge/Code%20Runner-orange?style=for-the-badge&logo=visual-studio-code&logoColor=white)  [![Visit My Page](https://img.shields.io/badge/Visit-My%20Page-orange?style=for-the-badge&logo=appveyor)](https://contact.andreyis.com)

This project contains various simulation schemes for stochastic processes, focusing on different **Test Functionals** and **Approximation Methods**.

Please refer to the page 58 ~66 of the paper [*High order polynomial regression approximation schemes in total variation for multidimensional diffusions*](https://www.overleaf.com/project/6801f9a43ca0501e11926ee2)

<!-- $\text{X\\_Y}$ \\_表示在github里面出现X_Y -->

---

## 📖 Table of Contents

*   **[Naming Convention](#naming-convention)**
    *   [C++ File Naming Convention](#c-file-naming-convention)
    *   [Data Source Naming Convention](#data-source-naming-convention)
*   **[Simulation Modules](#simulation-modules)**
    *   [1. 2DLT Series](#1-2dlt-series)
    *   [2. D1f Series](#2-d1f-series)
    *   [3. D1S Series (TV-distance Lower Bounds)](#3-d1s-series-tv-distance-lower-bounds)
    *   [4. D2f Series](#4-d2f-series)
    *   [5. D2S Series](#5-d2s-series)
    *   [6. Local Time (LT)](#6-local-time-lt)
    *   [7. Max Measure (M)](#7-max-measure-m)
    *   [8. Model 1 & Model 2 & Model 3(M1&M2&M3)](#8-model-1-model-2--model-3-m1-m2--m3)
* **[Calculation](#calculation)**

*   **[Development Environment](#development-environment)**
    *   [Compilation & Execution](#compilation--execution)
---
## Naming Convention

### C++ File Naming Convention
The file naming convention used throughout this project is as follows:

*   **`M1` / `M2`/ `M3` suffix** (e.g., `LTM1`, `D1fM2`,`2DD1SM3`):
    *   Indicates the specific **Model** used for simulation, as defined in the referenced paper.
    *   `M1` corresponds to [*Model 1*](#model-1).
    *   `M2` corresponds to [*Model 2*](#model-2).
    *   `M3` corresponds to [*Model 3*](#model-3).

*   **`_nobm` suffix** (e.g. `2DLTM3_nobm`):
    *   Indicates that the simulation is performed **without using a benchmark solution** ($\widetilde{X}^{\text{benchmark}}$).
    *   "bm" stands for benchmark. "nobm" means the exact solution or high-precision benchmark is not used/available for error calculation in this specific file.

*   **`D1` / `D2` prefix** (e.g., `D1fM1`, `D2S`):
    *   Indicates that the simulation is based on the functions $D_1(n, X)$ / $D_2(n, X)$.
    *   `D1` refers to the simulation of $D_1(n, X)$.
    *   `D2` refers to the simulation of $D_2(n, X)$.
    *   For specific formulas and definitions, please refer to the relevant sections in [Simulation Modules](#simulation-modules).

*   **`2D` prefix** (e.g., `2DLT`, `2DD1S`):
    *   Indicates that the simulation is performed in a **2-dimensional** setting.
    *   `2DD1S` refers to the simulation of $\mathrm{sgn}(D_1(n, X))$ in 2D, typically involving the **Sign** functional or related metrics.

 *   **`LT` prefix** (e.g., `LTM1`, `LTM2`):
        *   Indicates that the simulation focuses on the **Local Time** approximation.
        *   This prefix corresponds to the methods described in the [Local Time (LT)](#6-local-time-lt) section.
        *   Files with this prefix typically implement the discrete approximation of local time $L_n^z(X)$ and related test functionals.


*   **`M` prefix** (e.g., `MM1`, `MM2`):
    *   Indicates that the simulation focuses on the **Max Measure** (strong convergence).
    *   This prefix corresponds to the methods described in the [Max Measure (M)](#7-max-measure-m) section.
    *   Files with this prefix typically calculate the error defined by the maximum difference over the path, specifically, $\mathbb{E}[\max_{i=1,\ldots,n} \|X_{t_{i,n}} - \bar{X}^\theta_{t_{i,n}}\|^2]}$

*   **`_check` suffix** (e.g., `LTM1_check`):
    *   Represents **self-verification** of a function or scheme.  
    *   These files are used to confirm the behavior of a single scheme or function in isolation.  


*   **No suffix** (e.g., `LTM1`):
    *   Represents the **Error** between schemes. 
    *   These files typically calculate the difference between two methods using a specific formula to analyze **Convergence** or **Error Rates**.
    *   For Local Time simulations  `LTM1` and  `LTM2`, these files compute the error between the benchmark solution at the 1000th discretization point and the approximation scheme: $\mathbb{E}[F(\bar{X})] - \mathbb{E}[F(\bar{X}^\alpha)]$, where $\bar{X}$ is evaluated at $n=1000$ and $\bar{X}^\alpha$ represents the Euler, Milstein, or 1.5 order approximation.

*   **`_limit` suffix** (e.g., `D1fM1_limit`):
    *   Represents the **direct calculation of the limit**. 
    *   These files focus solely on computing the limiting value.  

*   **`_all` suffix** (e.g., `D1fM1_all`):
    *   Represents a **combined analysis**. 
    *   These files integrate both the error analysis (difference) and limit calculations into a single executable.
    *    :warning: *If a folder does not contain an `_all` file, then all simulations corresponding to `_limit` are included in the **no suffix** file.*

*  **`_var` suffix** (e.g., `LTM1_check_var`):
   *  It indicates files simulated using specific values. `LTM1_check_var` corresponds to `LTM1_check ` with $x_0=1, T=1, z=0.5, a=-1, b=1,\alpha=1.5$. They are mainly used to generate variance plots

*  **`_EFF` suffix** (e.g., `LTM1_EFF`):
   *  Represents the **Expected Functional Difference** calculation.
   *  These files compute $\mathbb{E}[F(\bar{X}) - F(\bar{X}^\alpha)]$, where $F$ is the test functional, $\bar{X}$ is the benchmark solution, and $\bar{X}^\alpha$ is the approximation scheme (Euler, Milstein, or 1.5 order).
   *  Used to evaluate the weak convergence error between the exact solution and the approximation.

### Data Source Naming Convention
The files in the `data_source/` directory follow a naming convention that links them to their generating C++ files and simulation parameters:

*   **`_100_1000_data` suffix**:
    *   Indicates that the dataset contains simulation results for a range of discretization points, specifically from $n=100$ to $n=1000$.

*   **File Prefix** (e.g., `D1fM1_check` in `D1fM1_check_100_1000_data.csv`):
    *   Corresponds directly to the name of the C++ file that generated the data (e.g., `D1fM1_check_100_1000_data.csv` is generated by `D1fM1_check.cpp`). 
    *   This convention applies to all variations including `_limit`, `_all`, `_var`, and `2D*` series (e.g., `2DD1SM1_100_1000_data.csv` from `2DD1SM1.cpp`), ensuring a clear mapping between the data and the code used to produce it.

---

<!-- \[
H_t^{112}(x,y)
= \frac{r_1^{2}\,r_2}{t^{3}\,s(x_2)^{4}\,s(x_1)^{2}}
-\frac{r_2}{t^{2}\,s(x_2)^{2}\,s(x_1)^{2}} 
\]

\[
H_t^{221}(x,y)
= \frac{r_2^{2}\,r_1}{t^{3}\,s(x_1)^{4}\,s(x_2)^{2}}
-\frac{r_1}{t^{2}\,s(x_1)^{2}\,s(x_2)^{2}} 
\] -->

## Simulation Modules

### 1. 2DLT Series
**Directory:** `2DLT/`

This module implements a test functional based on the **2D-Local Time** at a point $z$. It considers specific examples where:  

$$
\sigma(x) = \left( {s(x_2) \atop 0} \ \  {0 \atop s(x_1)} \right),\quad b(x) = \left( {s^{2}(x_2) \atop 0} \ \  {0 \atop s^{2}(x_1)} \right)
$$

Where,

$$
s(x) = 2 + \sin x, \quad a(x) = \left( {a_1(x) \atop a_2(x)}\right) = \left( {x_2 \atop -x_1}\right), \quad x_{\text{start}} = \left( {1 \atop 1}\right), \quad T = 1
$$

It implements various simulation schemes (Euler, Milstein, 1.5 order) for multidimensional processes.

Define discrete approximation of the local time at the point $z$ as

$$
L_n^z(X) = \sum_{k=1}^n \varphi_n(X_{t_{k,n}})h_n, \quad \varphi_n(x) = \sqrt{\frac{1}{2\pi h_n^\alpha}} e^{-\frac{(x-z)^2}{2h_n^\alpha}}, \quad \alpha \in (0, 2)
$$

We can use, for simplicity, a test functional based on the discrete approximation of the local time at the point $z$ for one component of the process

$$
F_n(\overline{X}) = f\left(L_n^z(\overline{X}^1)\right).\quad\cdots\text{(2DLT)}
$$

$$
f(x) = \arctan(x)
$$

*   **Simulation Plot:** [2D_plot](./2D/2D_plot.ipynb)
*   **Simulation Code:** [2D](./2D/)


### 2. D1f Series
**Directory:** `D1f/`

This module simulates the **D1f series** based on the mathematical definition of 

$$
D_1(n, X)=\sum_{k=1}^n \Delta^1(h_n, X_{t_{k-1}}, X_{t_k}) - \frac{1}{2} \sum_{k=1}^n \left[ \Delta^1(h_n, X_{t_{k-1}}, X_{t_k}) \right]^2
$$

Where,

$$
\begin{aligned}
\Delta^1(t, x, y) &= \frac{1}{4}b'(x)b(x)t^2 H_t^3(x, y) \\
&= \frac{1}{2}\sigma'(x)\sigma(x)^3 t^2 \left( \frac{(y-x-a(x)t)^3}{\sigma(x)^6 t^3} - 3\frac{(y-x-a(x)t)}{\sigma(x)^4 t^2} \right) \\
&= \frac{1}{2}\sigma'(x) \left( \frac{(y-x-a(x)t)^3}{\sigma(x)^3 t} - 3\frac{(y-x-a(x)t)}{\sigma(x)} \right)
\end{aligned}
$$

The functional is defined using a bounding function $f(x) = (x \wedge 0) \vee (-100)$ to ensure stability:  

$$
\bar{F}^E_n(X) = f(D_1(n, X))\quad\cdots\text{(D1f\\_check)}
$$

where $D_1(n, X)$ involves first-order increments.  

In this case the limit of $\mathbb{E} [\bar{F}_n^E(X^{[n]}) - \bar{F}_n^E(X^{1,[n]})]\quad\cdots\text{(D1f)}$ is: 

$$
\mathbb{E} \left[ f \left( I_T^0 + \frac{1}{2}\langle I^0 \rangle_T \right) \left( 1 - \exp \left( -I_T^0 - \frac{1}{2}\langle I^0 \rangle_T \right) \right) \right] \quad\cdots\text{(D1f\\_limit)}
$$

Where

$$
I_T^0 = \int_0^T \sqrt{\frac{3}{2}} |\sigma'|(X(s)) dW_s^{(3)}.  
$$

*   **Simulation Plot:** [D1f_plot](./D1f/D1f_plot.ipynb)
*   **Simulation Code:** [D1f](./D1f/)

### 3. D1S Series (TV-distance Lower Bounds)
**Directory:** `D1S/`, `2DD1S/`

This module focuses on test functionals for the lower bounds of the Total Variation (TV) distance.
It uses the sign function to define the functional:

$$
F_n^E(X) = \text{sgn}(D_1(n, X)) \quad\cdots\text{(D1S\\_check)}
$$

Where,

$$
\begin{aligned}
\Delta^1(t, x, y) &= \frac{1}{4}b'(x)b(x)t^2 H_t^3(x, y) \\
&= \frac{1}{2}\sigma'(x)\sigma(x)^3 t^2 \left( \frac{(y-x-a(x)t)^3}{\sigma(x)^6 t^3} - 3\frac{(y-x-a(x)t)}{\sigma(x)^4 t^2} \right) \\
&= \frac{1}{2}\sigma'(x) \left( \frac{(y-x-a(x)t)^3}{\sigma(x)^3 t} - 3\frac{(y-x-a(x)t)}{\sigma(x)} \right)
\end{aligned}
$$

This is used to show numerically that approximation rates for Euler schemes are not better than $h_n^0$.

The error is modeled as the expected difference between the benchmark $X^{\mathtt{benchmark}}$ and its order-0.5, order-1.0, and order-1.5 approximations, via the following difference-based simulation formula:  

$$
\begin{aligned}
\mathbb{E} \left[ F_n^E(X^{[n]}) - F_n^E(X^{0,[n]}) \right] &\approx \mathbb{E} F_n^E(X^{0,[n]}) (\mathcal{E}_n^1(X^{0,[n]}) - 1) \\
&\approx \mathbf{E} \mathrm{sgn} \left( \log \mathcal{E}_n^1(X^{0,[n]}) \right) (\mathcal{E}_n^1(X^{0,[n]}) - 1) = \mathbf{E} \left| \mathcal{E}_n^1(X^{0,[n]}) - 1 \right| \quad \dots \text{(D1S)} \\
&\approx d_{\mathrm{TV}} \left( \mathrm{Law}(X^{[n]}), \mathrm{Law}(X^{0,[n]}) \right),
\end{aligned}
$$

The limit can also be computed explicitly

$$
\mathbf{E} \left| \exp \left( -I_T^0 - \frac{1}{2}\langle I^0 \rangle_T \right) - 1 \right| \quad\cdots\text{(D1S\\_limit)}
$$

Where,

$$
I_T^0 = \int_0^T \sqrt{\frac{3}{2}} |\sigma'|(X(s)) dW_s^{(3)}.
$$

#### 2D Case Extension

In the 2-dimensional setting (Model 3), we extend the test functionals to:

**Test Functionals:**

$$
F_n^{E}(\overline{X})
=\textrm{sgn}\!\left(
\sum_{k=1}^{n}\Delta^{1}\!\left(h_n,X_{t_{k-1}},X_{t_k}\right)
-\frac{1}{2}\sum_{k=1}^{n}\left[\Delta^{1}\!\left(h_n,X_{t_{k-1}},X_{t_k}\right)\right]^2
\right)\quad\cdots (\texttt{2DD1S})
$$



**Definition of $\Delta_t^1$ (2D):**

$$
\Delta_t^{1}(x,y)
=\frac{1}{2} t^{2}s'(x_2)s(x_2)s^{2}(x_1)\,H_t^{112}(x,y)
+\frac{1}{2} t^{2}s'(x_1)s(x_1)s^{2}(x_2)\,H_t^{221}(x,y)
$$



*   **Simulation Plot:** [D1S_plot](./D1S/D1S_plot.ipynb), [2DD1S_plot](./2D/2D_plot.ipynb)
*   **Simulation Code:** [D1S](./D1S/), [2DD1S](./2D/)

### 4. D2f Series
**Directory:** `D2f/`

Define the $D_2$ function: 

$$
D_2(n, X) = \sum_{k=1}^n \Delta^2(h_n, X_{t_{k-1}}, X_{t_k})
$$

Where,

$$
\begin{aligned}
\Delta_t^2(x, y) &= \left( \frac{1}{2}a'(x)b(x) + \frac{1}{4}b'(x)a(x) + \frac{1}{8}b''(x)b(x) - \frac{1}{16}b'(x)^2 \right) t^2 H_t^2(x, y) \\
&\quad + \frac{1}{12}b''(x)b(x)^2 t^3 H_t^4(x, y) \\
&= \left( \frac{1}{2}a'(x)\sigma(x)^2 + \frac{1}{2}\sigma'(x)a(x)\sigma(x) + \frac{1}{4}\sigma''(x)\sigma(x)^3 \right) t^2 \left( \frac{(y-x-a(x)t)^2}{\sigma(x)^4 t^2} - \frac{1}{\sigma(x)^2 t} \right) \\
&\quad + \frac{1}{6} \left( \sigma''(x)\sigma(x)^5 + (\sigma'(x))^2\sigma(x)^4 \right) t^3 \left( \frac{(y-x-a(x)t)^4}{\sigma(x)^8 t^4} - 6\frac{(y-x-a(x)t)^2}{\sigma(x)^6 t^3} + \frac{3}{\sigma(x)^4 t^2} \right) \\
&= \left( \frac{1}{2}a'(x)\sigma(x) + \frac{1}{2}\sigma'(x)a(x) + \frac{1}{4}\sigma''(x)\sigma(x)^2 \right) \left( \frac{(y-x-a(x)t)^2}{\sigma(x)^3} - \frac{t}{\sigma(x)} \right) \\
&\quad + \frac{1}{6} \left( \sigma''(x)\sigma(x) + (\sigma'(x))^2 \right) \left( \frac{(y-x-a(x)t)^4}{\sigma(x)^4 t} - 6\frac{(y-x-a(x)t)^2}{\sigma(x)^2} + 3t \right).
\end{aligned}
$$

There is no theoretical prediction for the behavior of $\mathbb{E} \left[F_n^M(X^{[n]}) - F_n^M(X^{0,[n]})\right]$, the theoretical bound here $O(n^{1/2})$, but what we have observed in the earlier simulations is that the rate is very close to $O(n^{-1})$. 

In this set-up one may also consider

$$
\textrm{E} \left[\bar{F}_n^M(X^{[n]}) - \bar{F}_n^M(X^{1,[n]})\right] \quad\cdots \text{(D2f)}
$$

Where

$$
\bar{F}_n^M(X^{[n]}) := h_n^{-1/2} f(\Delta_2(n, X)) \quad\cdots\text{(D2f\\_check)}
$$

The limit is as before

$$
\mathbb{E} \left[ f \left(I_T^1\right) I_T^1 \right]. \quad\cdots \text{(D2f\\_limit)}
$$

In the case of the sign function this becomes $\mathbb{E}[|I_T^1|]$. 

*   **Simulation Plot:** [D2f_plot](./D2f/D2f_plot.ipynb)
*   **Simulation Code:** [D2f](./D2f/)

### 5. D2S Series
**Directory:** `D2S/`, `2DD2S/`

Similar to D1S, this module defines test functionals for lower bounds of TV-distance but for higher-order terms.

$$
F_n^M(X) = h_n^{-1/2}\text{sgn}(D_2(n, X)) \quad\cdots\text{(D2S\\_check)}
$$

where,

$$
\begin{aligned}
\Delta_t^2(x, y) &= \left( \frac{1}{2}a'(x)b(x) + \frac{1}{4}b'(x)a(x) + \frac{1}{8}b''(x)b(x) - \frac{1}{16}b'(x)^2 \right) t^2 H_t^2(x, y) \\
&\quad + \frac{1}{12}b''(x)b(x)^2 t^3 H_t^4(x, y) \\
&= \left( \frac{1}{2}a'(x)\sigma(x)^2 + \frac{1}{2}\sigma'(x)a(x)\sigma(x) + \frac{1}{4}\sigma''(x)\sigma(x)^3 \right) t^2 \left( \frac{(y-x-a(x)t)^2}{\sigma(x)^4 t^2} - \frac{1}{\sigma(x)^2 t} \right) \\
&\quad + \frac{1}{6} \left( \sigma''(x)\sigma(x)^5 + (\sigma'(x))^2\sigma(x)^4 \right) t^3 \left( \frac{(y-x-a(x)t)^4}{\sigma(x)^8 t^4} - 6\frac{(y-x-a(x)t)^2}{\sigma(x)^6 t^3} + \frac{3}{\sigma(x)^4 t^2} \right) \\
&= \left( \frac{1}{2}a'(x)\sigma(x) + \frac{1}{2}\sigma'(x)a(x) + \frac{1}{4}\sigma''(x)\sigma(x)^2 \right) \left( \frac{(y-x-a(x)t)^2}{\sigma(x)^3} - \frac{t}{\sigma(x)} \right) \\
&\quad + \frac{1}{6} \left( \sigma''(x)\sigma(x) + (\sigma'(x))^2 \right) \left( \frac{(y-x-a(x)t)^4}{\sigma(x)^4 t} - 6\frac{(y-x-a(x)t)^2}{\sigma(x)^2} + 3t \right).
\end{aligned}
$$

Similarly, we expect to observe that

$$
\begin{aligned}
&\mathbb{E} \left[ F_n^M(X^{[n]}) - F_n^M(X^{1,[n]}) \right] \quad\cdots\text{(D2S)} \\ 
&\approx \mathbf{E} \text{sgn} \left( \log \mathcal{E}_n^2(X^{0,[n]}) \right) (\mathcal{E}_n^2(X^{0,[n]}) - 1) \\
&\to \mathbf{E} \text{sgn}(\log(1 + I_T^1))(I_T^1) \quad\cdots\text{(D2S\\_limit)}
\end{aligned}
$$

Where,

$$
I_T^1 = \sqrt{2} \int_0^T |c_2^2||b|^{-1}(X(s)) dW_s^{(2)} + \sqrt{4! } \int_0^T |c_4^2| b^{-2}(X(s)) dW_s^{(4)}
$$

$$
\begin{aligned}
c_4^2(x) &= \frac{1}{12}b''(x)b(x)^2 \\
c_2^2(x) &= \frac{1}{2}a'(x)b(x) + \frac{1}{4}b'(x)a(x) + \frac{1}{8}b''(x)b(x) - \frac{1}{16}b'(x)^2. 
\end{aligned}
$$

Furthermore,

$$
\mathbb{E} \left[F_n^M(X^{[n]}) - F_n^M(X^{2,[n]})\right] = O(n^{-1/2}).
$$

It is used to analyze the Milstein scheme's approximation rates. 

#### 2D Case Extension

In the 2-dimensional setting (Model 3), we extend the test functional for higher-order terms:

**Test Functional:**

$$
F_n^{M}(\overline{X})
=h_n^{-1/2}\textrm{sgn}\!\left(
\sum_{k=1}^{n}\Delta^{2}\!\left(h_n,X_{t_{k-1}},X_{t_k}\right)
\right)\quad\cdots (\texttt{2DD2S})
$$

**Definition of $\Delta_t^2$ (2D):**

<!-- $$
\begin{aligned}
\Delta_t^{2}(x,y)
&=\frac{t^{2}}{2}\,\partial_{1}a_{1}(x)s^2(x_{2})\,H_t^{11}(x,y)+\frac{t^{2}}{2}\,\partial_{2}a_{1}(x)s^2(x_{1})\,H_t^{12}(x,y) \\
&\quad+\frac{t^{2}}{2}\,\partial_{1}a_{2}(x)s^2(x_{2})\,H_t^{12}(x,y)+\frac{t^{2}}{2}\,\partial_{2}a_{2}(x)s^2(x_{1})\,H_t^{22}(x,y) \\
&\quad+\frac{t^{2}}{4}\,s^2(x_{2})s'(x_{2})a_{2}(x)\,H_t^{11}(x,y)+\frac{t^{2}}{4}\,s^2(x_{1})s'(x_{1})a_{1}(x)\,H_t^{22}(x,y) \\
&\quad+\frac{t^{2}}{4}\,s''(x_{2})s(x_{2})s(x_{1})^{2}\,H_t^{11}(x,y)+\frac{t^{2}}{4}\,s''(x_{1})s(x_{1})s(x_{2})^{2}\,H_t^{22}(x,y) \\
&\quad+\frac{t^{2}}{8}\,[s'(x_{2})s(x_{1})]^{2}\,H_t^{11}(x,y)+\frac{t^{2}}{8}\,[s'(x_{1})s(x_{2})]^{2}\,H_t^{22}(x,y) \\
&\quad-\frac{t^{2}}{4}\,s'(x_{1})s(x_{1})s'(x_{2})s(x_{2})\,H_t^{12}(x,y) \\
&\quad+\frac{t^{3}}{24}\,(s'(x_{2}))^{2}s(x_{1})^{2}s(x_{2})^{2}\,H_t^{1111}(x,y)+\frac{t^{3}}{24}\,(s'(x_{1}))^{2}s(x_{1})^{2}s(x_{2})^{2}\,H_t^{2222}(x,y) \\
&\quad+\frac{t^{3}}{12}\,s'(x_{1})s'(x_{2})s(x_{1})s(x_{2})^{3}\,H_t^{1112}(x,y)+\frac{t^{3}}{12}\,s'(x_{1})s'(x_{2})s(x_{1})^{3}s(x_{2})\,H_t^{1222}(x,y) \\
&\quad+\left[\frac{t^{3}}{6}\,s''(x_{1})s(x_{1})s(x_{2})^{4}+\frac{t^{3}}{6}\,s''(x_{2})s(x_{2})s(x_{1})^{4} \right. \\
&\quad\quad\left. +\frac{t^{3}}{24}\,(s'(x_{2}))^{2}s(x_{1})^{4}+\frac{t^{3}}{24}\,(s'(x_{1}))^{2}s(x_{2})^{4}\right] H_t^{1122}(x,y).
\end{aligned}
$$ -->

$$
\begin{aligned}
\Delta_t^{2}(x,y)
&= \frac{t^{2}}{2} \partial_{1}a_{1}(x)s^2(x_{2}) H_t^{11}(x,y)+ \frac{t^{2}}{2} \partial_{2}a_{1}(x)s^2(x_{1}) H_t^{12}(x,y) \\
&\quad + \frac{t^{2}}{2} \partial_{1}a_{2}(x)s^2(x_{2}) H_t^{12}(x,y)+ \frac{t^{2}}{2} \partial_{2}a_{2}(x)s^2(x_{1}) H_t^{22}(x,y)
\end{aligned}
$$
$$
\begin{aligned}
&\quad + \frac{t^{2}}{4} s^2(x_{2})s'(x_{2})a_{2}(x) H_t^{11}(x,y)+ \frac{t^{2}}{4} s^2(x_{1})s'(x_{1})a_{1}(x) H_t^{22}(x,y) \\
&\quad + \frac{t^{2}}{4} s''(x_{2})s(x_{2})s(x_{1})^{2} H_t^{11}(x,y)+ \frac{t^{2}}{4} s''(x_{1})s(x_{1})s(x_{2})^{2} H_t^{22}(x,y)
\end{aligned}
$$
$$
\begin{aligned}
&\quad + \frac{t^{2}}{8} [s'(x_{2})s(x_{1})]^{2} H_t^{11}(x,y)+ \frac{t^{2}}{8} [s'(x_{1})s(x_{2})]^{2} H_t^{22}(x,y) \\
&\quad - \frac{t^{2}}{4} s'(x_{1})s(x_{1})s'(x_{2})s(x_{2}) H_t^{12}(x,y)
\end{aligned}
$$
$$
\begin{aligned}
&\quad + \frac{t^{3}}{24} (s'(x_{2}))^{2}s(x_{1})^{2}s(x_{2})^{2} H_t^{1111}(x,y)+ \frac{t^{3}}{24} (s'(x_{1}))^{2}s(x_{1})^{2}s(x_{2})^{2} H_t^{2222}(x,y)
\end{aligned}
$$
$$
\begin{aligned}
&\quad + \frac{t^{3}}{12} s'(x_{1})s'(x_{2})s(x_{1})s(x_{2})^{3} H_t^{1112}(x,y)+ \frac{t^{3}}{12} s'(x_{1})s'(x_{2})s(x_{1})^{3}s(x_{2}) H_t^{1222}(x,y)
\end{aligned}
$$
$$
\begin{aligned}
&\quad + \frac{t^{3}}{6} s''(x_{1})s(x_{1})s(x_{2})^{4} H_t^{1122}(x,y)+ \frac{t^{3}}{6} s''(x_{2})s(x_{2})s(x_{1})^{4} H_t^{1122}(x,y) \\
&\quad + \frac{t^{3}}{24} (s'(x_{2}))^{2}s(x_{1})^{4} H_t^{1122}(x,y)+ \frac{t^{3}}{24} (s'(x_{1}))^{2}s(x_{2})^{4} H_t^{1122}(x,y)
\end{aligned}
$$


*   **Simulation Plot:** [D2S_plot](./D2S/D2S_plot.ipynb), [2DD2S_plot](./2D/2D_plot.ipynb)
*   **Simulation Code:** [D2S](./D2S/), [2DD2S](./2D/)

### 6. Local Time (LT)
**Directory:** `LT/`

This module implements a test functional based on the **Local Time** at a point $z$. 

Define discrete approximation of the local time at the point $z$ as

$$
L_n^z(X) = \sum_{k=1}^n \varphi_n(X_{t_{k,n}})h_n \quad\cdots\text{(LT\\_check)}, \quad \varphi_n(x) = \sqrt{\frac{1}{2\pi h_n^\alpha}} e^{-\frac{(x-z)^2}{2h_n^\alpha}}, \quad \alpha \in (0, 2)
$$

The test functional used is $F_n(\overline{X}^\theta) = \arctan(L_n^z(\overline{X}^\theta))$.

Particular choice:  $f(x) = \arctan(x)$. We will compare

$$
\mathbb{E}[F_n(\overline{X})] - \mathbb{E}[F_n(\overline{X}^\alpha)] \quad\cdots\text{(LT)}
$$

*   **Simulation Plot:** [LT_plot](./LT/LT_plot.ipynb)
*   **Simulation Code:** [LT](./LT/)

### 7. Max Measure (M)
**Directory:** `Max/`

This module measures errors using the **maximum measure** (strong convergence).
It compares different schemes (Euler, Milstein, 1.5) using:  

$$
\mathbb{E} \left[ \max_{i=1,\ldots,n} \|X_{t_{i,n}} - \bar{X}^\theta_{t_{i,n}}\|^2 \right] \quad\cdots\text{(M)}
$$

*   **Simulation Plot:** [MM_plot](./Max/MM_plot.ipynb)
*   **Simulation Code:** [Max](./Max/)

### 8. Model 1, Model 2 & Model 3 (M1, M2 & M3)
**Directory:** `*M1/M2/M3`

#### Model 1

*  **M1 benchmark**:

$$
    X_t^{\mathtt{benchmark}} = \sinh \left( \ln(x + \sqrt{x^2 + 1}) + \frac{ab}{2}t + b W_t \right) 
$$

* **The SDE for $X_t$ of M1**:

$$
dX_t = \left( \frac{b^2}{2}X_t + \frac{ab}{2}\sqrt{X_t^2+1} \right) dt + b \sqrt{X_t^2+1} \, dW_t,
$$

Where,

$$
a(X)=\frac{b^2}{2}X_t + \frac{ab}{2}\sqrt {X_t^2+1}\ ,\sigma(X)= b \sqrt{X_t^2+1}
$$

#### Model 2

*  **M2 benchmark**:

$$
Y_t = e^{at}Y_0 + b \int_0^t e^{a(t-s)} \, dW_s = e^{at} \mathrm{asinh(X_0)}+ b \int_0^t e^{a(t-s)} \, dW_s.
$$

And then,

$$
    X_t^{\mathtt{benchmark}} = \sinh(Y_t)
$$

*  **The SDE for $X_t$ of M2**:

$$
dX_t = \left[ a \sqrt{1 + X_t^2} \mathrm{asinh}(X_t) + \frac{b^2}{2} X_t \right] dt + b \sqrt{1 + X_t^2} \, dW_t.
$$

where,

$$
a(X)=a \sqrt{1 + X_t^2} \mathrm{asinh}(X_t) + \frac{b^2}{2} X_t\, ,\sigma(X) = b \sqrt{1 + X_t^2}
$$

#### Model 3

*   **System Setup (Multivariate Case)**:

In the multivariate case, we consider specific example

$$
\sigma(x) = \begin{pmatrix} s(x_2) & 0 \\ 0 & s(x_1) \end{pmatrix}, \quad b(x) = \begin{pmatrix} s^2(x_2) & 0 \\ 0 & s^2(x_1) \end{pmatrix}
$$

<!-- $\color{red}{\text{Explicit Definition of } a \text{ and } s}$ -->

$$
s(x) = 2 + \sin x, \quad a(x) = \begin{pmatrix} a_1(x) \\ a_2(x) \end{pmatrix} = \begin{pmatrix} x_2 \\ -x_1 \end{pmatrix}, \quad x_{start} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}, \quad T = 1
$$
*   **Simulation schemes**:


$w_1, w_2 \sim N(0, t)$ are independent

$$
\begin{aligned}
\Delta_t^{2}(x,y)
&= \frac{t^{2}}{2} \partial_{1}a_{1}(x)s^2(x_{2}) H_t^{11}(x,y)+ \frac{t^{2}}{2} \partial_{2}a_{1}(x)s^2(x_{1}) H_t^{12}(x,y) \\
&\quad + \frac{t^{2}}{2} \partial_{1}a_{2}(x)s^2(x_{2}) H_t^{12}(x,y)+ \frac{t^{2}}{2} \partial_{2}a_{2}(x)s^2(x_{1}) H_t^{22}(x,y)
\end{aligned}
$$
$$
\begin{aligned}
&\quad + \frac{t^{2}}{4} s^2(x_{2})s'(x_{2})a_{2}(x) H_t^{11}(x,y)+ \frac{t^{2}}{4} s^2(x_{1})s'(x_{1})a_{1}(x) H_t^{22}(x,y) \\
&\quad + \frac{t^{2}}{4} s''(x_{2})s(x_{2})s(x_{1})^{2} H_t^{11}(x,y)+ \frac{t^{2}}{4} s''(x_{1})s(x_{1})s(x_{2})^{2} H_t^{22}(x,y)
\end{aligned}
$$
$$
\begin{aligned}
&\quad + \frac{t^{2}}{8} [s'(x_{2})s(x_{1})]^{2} H_t^{11}(x,y)+ \frac{t^{2}}{8} [s'(x_{1})s(x_{2})]^{2} H_t^{22}(x,y) \\
&\quad - \frac{t^{2}}{4} s'(x_{1})s(x_{1})s'(x_{2})s(x_{2}) H_t^{12}(x,y)
\end{aligned}
$$
$$
\begin{aligned}
&\quad + \frac{t^{3}}{24} (s'(x_{2}))^{2}s(x_{1})^{2}s(x_{2})^{2} H_t^{1111}(x,y)+ \frac{t^{3}}{24} (s'(x_{1}))^{2}s(x_{1})^{2}s(x_{2})^{2} H_t^{2222}(x,y)
\end{aligned}
$$
$$
\begin{aligned}
&\quad + \frac{t^{3}}{12} s'(x_{1})s'(x_{2})s(x_{1})s(x_{2})^{3} H_t^{1112}(x,y)+ \frac{t^{3}}{12} s'(x_{1})s'(x_{2})s(x_{1})^{3}s(x_{2}) H_t^{1222}(x,y)
\end{aligned}
$$
$$
\begin{aligned}
&\quad + \frac{t^{3}}{6} s''(x_{1})s(x_{1})s(x_{2})^{4} H_t^{1122}(x,y)+ \frac{t^{3}}{6} s''(x_{2})s(x_{2})s(x_{1})^{4} H_t^{1122}(x,y) \\
&\quad + \frac{t^{3}}{24} (s'(x_{2}))^{2}s(x_{1})^{4} H_t^{1122}(x,y)+ \frac{t^{3}}{24} (s'(x_{1}))^{2}s(x_{2})^{4} H_t^{1122}(x,y)
\end{aligned}
$$
where

$$
w_t^{111} = (w_1)^3 - 3tw_1,\quad
w_t^{222} = (w_2)^3 - 3tw_2,\quad
w_t^{112} = \big((w_1)^2 - t\big)w_2,\quad
w_t^{122} = \big((w_2)^2 - t\big)w_1
$$

Please refer to the page 59 ~ 61 of the [*paper*](https://www.overleaf.com/project/6801f9a43ca0501e11926ee2)

---

## Calculation
**Directory:** `Calculation/`

This module contains the detailed mathematical derivations for the higher-order Hermite polynomials used in the simulations, specifically $H^{112}_t(x,y)$, $H^{221}_t(x,y)$ of $\Delta_t^1$ and $H^{1222}_t(x,y)$ of $\Delta_t^2$ and so on.

*   **Compile File `.tex` :** `higher_order_hermite_calculation.tex`
*  **PDF File :** *[higher_order_hermite_calculation](./Calculation/out/higher_order_hermite_calculation.pdf)*
*   **Content:**
    *   Setup and definitions of $\sigma(x)$ and $b(x)$.
    *   Derivation of centralized variables and standardized variables.
    *   Step-by-step calculation of $H^{i,j,k, \cdots}_t(x,y)$ of $\Delta_t^1$ and $\Delta_t^2$ based on the general formula.(please refer to the page 53, Eq. A.5 or A.7 of the *[PAPER](https://www.overleaf.com/project/6801f9a43ca0501e11926ee2)*


---

## Development Environment

This project is developed using **Visual Studio Code** IDE.

### Compilation & Execution
To compile and run the C++ simulations, we recommend using the **Code Runner** extension.

*   **OpenMP Support**: The C++ code utilizes **OpenMP** for parallel processing to accelerate simulations.
*   **Configuration**:
    *   Essential configuration settings, including include paths for **Eigen** and **OpenMP**, are defined in `.vscode/c_cpp_properties.json` and `.vscode/settings.json`.
    *   These settings are **mandatory** for successful compilation, as they specify the necessary compiler flags and library paths.

Please refer to the `.vscode/` directory for specific include paths and environment settings.
