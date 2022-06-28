# Homework 4

**周雨豪  2018013399  软件92**



## 1 采样方法

**(a) 证明拒绝采样得到的分布就是原分布**

由拒绝采样定义得知存在 $k$ 满足 $k\cdot q(x)\geq p(x)$

考虑：$x\sim q$，$u\sim U[0,k\cdot q(x)]$，如果 $u\gt p(x)$ 则拒绝 $x$

有 $p(x|\text{accept})\propto q(x)p(\text{accept}|x)=q(x)(\frac{p(x)}{k\cdot q(x)})=\frac{p(x)}{k}$

因此拒绝采样得到的分布就是原分布

**(b) 证明 Gibbs 采样法的过程满足细致平衡**

由 Gibbs 采样法得知 $x_{\lnot k}^*=x_{\lnot k}$，$p(x)=p(x_{\lnot k})p(x_k|x_{\lnot k})$

$\frac{\pi_jQ_{ji}}{\pi_iQ_{ij}}=\frac{p(x^*)q(x|x^*)}{p(x)q(x^*|x)}=\frac{p(x_{\lnot k}^*)p(x_k^*|x_{\lnot k}^*)p(x_k|x_{\lnot k}^*)}{p(x_{\lnot k})p(x_k|x_{\lnot k})p(x_k^*|x_{\lnot k})}=1$

证得 $\pi_jQ_{ji}=\pi_iQ_{ij}$，因此 Gibbs 采样法的过程满足细致平衡



## 2 Baum-Welch 算法

**(a)**



**(b)**

**(c)**

**(d)**

**(e)**