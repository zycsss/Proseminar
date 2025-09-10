# Important Functions in the Paper

## 1. **Compression Complexity Function**  
**Equation**:  
$$\eta(\beta_k, \epsilon) = e^{\beta_k \cdot \epsilon} - e^{\epsilon}$$  
**Explanation**:  
Models the computational complexity of lossless compression at sensor \(k\). Complexity increases exponentially with the compression ratio \($\beta_k$\), where \($\epsilon$\) is a compression algorithm-specific constant.  

---

## 2. **Compression Time**  
**Equation**:  
$$T_k^{\text{comp}}(\beta_k) = \frac{D_k \cdot \eta(\beta_k, \epsilon)}{f_k^{\text{S}}}$$  
**Explanation**:  
Time for sensor \($k$\) to compress \($D_k$\) bits of data. Depends on compression complexity \($\eta(\beta_k, \epsilon)$\), data volume \($D_k$\), and the sensor's compression speed \($f_k^{\text{S}}$\).  

---

## 3. **Achievable Data Rate**  
**Equation**:  
$$r_k(b_k) = b_k \cdot \log_2\left(1 + \frac{|H_k|^2 \cdot p_k}{N_0 \cdot b_k}\right)$$  
**Explanation**:  
Shannon capacity formula for OFDMA transmission. Represents the data rate of sensor \($k$\) given bandwidth \($b_k$\), transmit power \($p_k$\), channel gain \($H_k$\), and noise spectral density \($N_0$\).  

---

## 4. **Transmission Time**  
**Equation**:  
$$T_k^{\text{tr}}(\beta_k, b_k) = \frac{D_k}{\beta_k \cdot r_k(b_k)}$$  
**Explanation**:  
Time for sensor \($k$\) to transmit compressed data. Inversely proportional to the compression ratio \($\beta_k$\) and achievable data rate \($r_k(b_k)$\).  

---

## 5. **DT Processing Time**  
**Equation**:  
$$T_k^{\text{DT}}(f_k^{\text{DT}}) = \frac{D_k \cdot c_k}{f_k^{\text{DT}}}$$  
**Explanation**:  
Time for the cloud server to process data from sensor \($k$\). Depends on computational complexity \($c_k$\) and allocated processing speed ($f_k^{\text{DT}}$\).  

---

## 6. **Total DT Synchronization Time per Sensor**  
**Equation**:  
$$
T_k^{\text{total}}(\beta_k, b_k, f_k^{\text{DT}}) = T_k^{\text{comp}} + T_k^{\text{tr}} + T_k^{\text{DT}}
$$


**Explanation**:  
Total time for sensor \(k\) to compress, transmit, and process data. The global DT synchronization time is the maximum \($T_k^{\text{total}}$\) across all sensors.  

---

## 7. **Follower’s Cost Function**  
**Equation**:  
$$U_k^{\text{f}}(\beta_k; b_k) = T_k^{\text{comp}}(\beta_k) + T_k^{\text{tr}}(\beta_k, b_k)$$  
**Explanation**:  
Each sensor \(k\) minimizes this cost to balance compression time and transmission time for a given bandwidth \($b_k$\).  

---

## 8. **Leader’s Cost Function**  
**Equation**:  
$$
U^{\text{l}}(\mathbf{b}, \mathbf{f}; \boldsymbol{\beta}) = \max_{k \in \mathcal{K}} T_k^{\text{total}}(\beta_k, b_k, f_k^{\text{DT}})
$$
**Explanation**:  
The BS minimizes the worst-case DT synchronization time by optimizing bandwidth allocation \($\mathbf{b}$\) and processing speed \($\mathbf{f}$\).  

---

## 9. **Reformulated Follower’s Objective**  
**Equation**:  
$$f_k(\mathbf{x}_k; b_k) = \frac{D_k}{f_k^{\text{S}}} \cdot e^{\beta_k \cdot \epsilon} + \frac{D_k}{r_k(b_k)} \cdot \mu_k$$  
**Explanation**:  
Convex reformulation using auxiliary variable \($\mu_k = 1/\beta_k$\) to derive closed-form solutions. Combines compression and transmission costs.  

---

## 10. **Optimal Compression Ratio (\(\beta_k^*\)) and Auxiliary Variable (\(\mu_k^*\))**  
**Equations**:  
$$\beta_k^* = \min\left\{\frac{1}{\epsilon} \ln\left(\frac{f_k^{\text{S}}}{\epsilon D_k} (\lambda_{1,k} + \lambda_{2,k})\right), \beta^{\max}\right\}$$  
$$\mu_k^* = \sqrt{\lambda_{2,k} \cdot r_k(b_k) / D_k}$$  
**Explanation**:  
Closed-form solutions for sensor \(k\)’s best response. Derived via Lagrangian duality, where \($\lambda_{1,k}$\) and \($\lambda_{2,k}$\) are dual variables for constraints.  

---

## 11. **Lagrangian Function**  
**Equation**:  
$$L_k(\boldsymbol{\lambda}_k, \mathbf{x}_k; b_k) = f_k(\mathbf{x}_k; b_k) + \boldsymbol{\lambda}_k^{\intercal} \mathbf{g}(\mathbf{x}_k)$$  
**Explanation**:  
Combines the follower’s objective with dual variables ($\boldsymbol{\lambda_k}$) and constraints \($\mathbf{g}(\mathbf{x}_k)$\). Used to solve the dual problem for equilibrium.  

---

## 12. **Sub-Gradient Update Rule**  
**Equation**:  
$$\boldsymbol{\lambda}_k^{(l+1)} = \left[\boldsymbol{\lambda}_k^{(l)} + \mathbf{h}_k^{(l)} \circ \mathbf{g}\left(\mathbf{x}_k^*(\boldsymbol{\lambda}_k^{(l)})\right)\right]^+$$  
**Explanation**:  
Iterative update for dual variables using adaptive step size \($\mathbf{h}_k^{(l)}$\). Ensures convergence to the Stackelberg equilibrium.  

---

## Interesting CVX functions
- Leader problem:
    - Relative Entropy function: `rel_entr(x,y) = x .* log(x/y)`
    - Inverse of a product: `prop_inv(x,y)`
    - Positive inverse: `inv_pos(x)`
- Reformulation
    $$T_k^t >= D_k /(\beta \cdot r(b))$
- $T_k^t >= \frac{D_k} {\beta}\times\frac{1}{b \cdot log(1+c/b)}$
    $b \cdot \log(1+c/b) >= \frac{D_k}{\beta \times T_k^t }$
- Relative entropy function `rel_entr(x,y)`
    - `-rel_entr(x,y) = - (x .* log(x/y)) = x .* log(y/x)`
        - `x = b`, ,`y = b+c`