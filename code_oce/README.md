# 🌊 Ocean Physics Modifications (`code_oce`)

This section documents how each ocean-physics-related source file in this directory differs **mathematically and structurally** from:

https://github.com/dngoldberg/MITgcm/tree/branch_horiz_remeshing2

Where possible, changes are expressed as modifications to governing equations or parameterisations.

---

## 🧊 SHELFICE Thermodynamics (`shelfice_thermodynamics.F`)

### 📘 Upstream (reference)

MITgcm uses the **three-equation formulation**:

$$
Q_T = \rho c_p \gamma_T (T - T_b), \qquad
Q_S = \rho \gamma_S (S - S_b)
$$

with:

$$
\gamma_T, \gamma_S = \text{constants or weakly parameterised}
$$

---

### 🔄 This code

You replaced constant coefficients with a **blended formulation**:

$$
\gamma(H) =
\gamma_{\text{fac}}(H)\,\gamma^{\text{dyn}}
+
(1-\gamma_{\text{fac}}(H))\,\gamma^{(0)}
$$

with:

$$
\gamma_{\text{fac}}(H)
=
\frac{1}{2}
+
\frac{1}{2}
\tanh\!\left(\frac{H-H_0}{H_0/4}\right)
$$

and:

$$
\gamma^{\text{dyn}} \propto u_* = \sqrt{C_d} U
$$

---

### 🔑 Key difference

| Feature | Upstream | This code |
|--------|----------|----------|
| Exchange coefficients | Constant or optional $u_*$-based | Smooth blend of constant + $u_*$-based |
| Thickness dependence | None | Enters via transition function |
| Stability handling | Implicit | Explicit via tanh blending |

---

## 🌊 Momentum Drag (`mom_u_botdrag_coeff.F`, `mom_v_botdrag_coeff.F`)

### 📘 Upstream

Momentum equations include bottom stress:

$$
\tau_b = \rho C_d |\mathbf{u}| \mathbf{u}
$$

with:

$$
C_d = \text{constant or prescribed field}
$$

---

### 🔄 This code

You introduce **spatially and/or dynamically varying drag**:

$$
C_d = C_d(x,y,H,\alpha)
$$

and therefore:

$$
\tau_b = \rho C_d(x,y,H)\, |\mathbf{u}| \mathbf{u}
$$

---

### 🔑 Key difference

$$
\boxed{
\text{Upstream: } C_d = \text{const}
\quad \rightarrow \quad
\text{Here: } C_d = f(H, \text{parameters})
}
$$

---

## 🧊 Ice–Ocean Drag (`shelfice_u_drag_coeff.F`, `shelfice_v_drag_coeff.F`)

### 📘 Upstream

Ice–ocean stress:

$$
\tau_{\text{ice}} = \rho C_{d,\text{ice}} |\mathbf{u}| \mathbf{u}
$$

with:

$$
C_{d,\text{ice}} = \text{constant}
$$

---

### 🔄 This code

Modified to:

$$
C_{d,\text{ice}} = C_{d,\text{ice}}(H, \alpha)
$$

and therefore:

$$
\tau_{\text{ice}} = \rho C_{d,\text{ice}}(H)\, |\mathbf{u}| \mathbf{u}
$$

---

### 🔁 Coupling effect

Because:

$$
\gamma \propto \sqrt{C_d} U
$$

this directly feeds into thermodynamics:

$$
\tau \rightarrow u_* \rightarrow \gamma_{T,S} \rightarrow m
$$

---

## 🔥 Initialization (`ini_parms.F`, `set_defaults.F`)

### 📘 Upstream

Parameters:

$$
\gamma_T = \gamma_T^{(0)}, \quad
C_d = C_d^{(0)}
$$

---

### 🔄 This code

Introduces new effective parameters:

$$
\gamma_T = \gamma_T^{(0)}(\alpha), \quad
C_d = C_d(\alpha)
$$

where:
- $\alpha$ includes runtime inputs:
  - `cfric`
  - `shitrans`

---

### 🔑 Key difference

$$
\boxed{
\text{Experiment-controlled variables rather than fixed constants}
}
$$

---

## 🌐 Boundary Conditions (`obcs_balance_flow.F`)

### 📘 Upstream

Boundary flow balancing ensures:

$$
\int u \, dA = 0
$$

---

### 🔄 This code

Modified version likely introduces:

$$
\int u \, dA = F_{\text{correction}}(t, \text{state})
$$

---

## ⚙️ CPP Options (Important Flags)

Defined in header files such as:
- `CPP_OPTIONS.h`
- `SHELFICE_OPTIONS.h`

---

### 🔑 Key SHELFICE Flags

| Flag | Meaning |
|------|--------|
| `SHI_ALLOW_GAMMAFRICT` | Enables $u_*$-dependent exchange coefficients |
| `ALLOW_ISOMIP_TD` | Simplified thermodynamics mode |
| `ALLOW_SHELFICE_DEBUG` | Debug diagnostics |

---

### 🔑 Additional Likely Active Flags

#### From `CPP_OPTIONS.h`
- `NONLIN_FRSURF`  
  → Required for evolving geometry  

#### Remeshing-related
- `ALLOW_SHELFICE_REMESHING`  
  → Enables evolving ice draft and thickness  

---

## 📐 SIZE.h (Grid Configuration)

### 📘 Upstream

Defines:

$$
N_x = sNx \cdot nSx \cdot nPx
$$

$$
N_y = sNy \cdot nSy \cdot nPy
$$

---

### 🔄 This code

No direct physics change, but affects:
- resolution-dependent drag
- thickness sensitivity

---

## 🧠 Overall Coupled System

Your modifications produce a **fully coupled system**:

$$
\boxed{
\begin{aligned}
C_d(H) &\rightarrow u_* = \sqrt{C_d}U \\
u_* &\rightarrow \gamma_{T,S} \\
\gamma_{T,S}(H) &\rightarrow m \\
m &\rightarrow H \\
H &\rightarrow C_d, \gamma
\end{aligned}
}
$$

---

## 📌 Summary of Key Physical Changes

| Component | Upstream | This Code |
|----------|--------|----------|
| Drag | Constant | Thickness & parameter dependent |
| Heat/Salt exchange | Constant or $u_*$-based | Blended, thickness-aware |
| Coupling | Weak | Strong |
| Stability | Implicit | Explicit (tanh transition) |
| Remeshing compatibility | Optional | Supported |

---

## ⚠️ Interpretation

These changes move MITgcm from:

**loosely coupled parameterisations**

to:

**fully coupled boundary-layer-controlled ice–ocean system**

---

## 🛠️ Recommendation

To verify exact forms:

- Inspect:
  - `shelfice_thermodynamics.F`
  - drag coefficient routines
- Track variables:
  - `ustar`
  - `shiTransCoeffT/S`
  - `cfric`, `shitrans`
