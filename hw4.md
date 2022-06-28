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
$$
\begin{align*}
&J(\theta)=\mathcal{L}(q,\theta)=\sum\limits_Zq(Z)\log(\frac{P(X,Z|\theta)}{q(Z)}) \\
&=\sum\limits_Zq(Z)\log(P(X,Z|\theta))-\sum\limits_Z q(Z)\log(q(Z)) \\
&\Rightarrow \mathop{\arg\max}\limits_\theta J(\theta)=\mathop{\arg\max}\limits_\theta\sum\limits_Z q(Z)\log(P(X,Z|\theta))\\
&又\ q(Z)=P(Z|X,\theta^{(i)})=\frac{P(X,Z|\theta^{(i)})}{P(X|\theta^{(i)})}，P(X|\theta^{(i)})与\theta 无关\\
&因此\ \mathop{\arg\max}\limits_\theta J(\theta)=\mathop{\arg\max}\limits_\theta\sum\limits_ZP(X,Z|\theta^{(i)})\log P(X,Z|\theta)
\end{align*}
$$
**(b)**
$$
\begin{align*}
P(X,Z|\theta)&=\pi_{z_1}\prod\limits_{t=1}\limits^{T-1}a_{z_{t}z_{t+1}}\prod\limits_{t=1}\limits^{T}b_{z_t}(x_t) \\
Q(\theta,\theta^{(i)})&=\sum\limits_{Z}P(X,Z|\theta^{(i)})\log P(X,Z|\theta)\\
&=\sum\limits_{Z}P(X,Z|\theta^{(i)})\log \pi_{z_1}\prod\limits_{t=1}\limits^{T-1}a_{z_{t}z_{t+1}}\prod\limits_{t=1}\limits^{T}b_{z_t}(x_t) \\
&=\sum\limits_{Z}P(X,Z|\theta^{(i)})\log \pi_{z_1}+
\sum\limits_{Z}P(X,Z|\theta^{(i)})\sum\limits_{t=1}\limits^{T-1}\log a_{z_{t}z_{t+1}}+
\sum\limits_{Z}P(X,Z|\theta^{(i)})\sum\limits_{t=1}\limits^{T}\log b_{z_t}(x_t)
\end{align*}
$$
**(c)**
$$
\begin{align*}
\sum\limits_Z P(X,Z|\theta^{(i)})\log\pi_{z_1}&=\sum\limits_{z_1}\dots\sum\limits_{z_T}P(X,z_1,\dots,z_T|\theta^{(i)})\\
&=\sum\limits_{z_1}\log\pi_{z_1}P(X,z_1|\theta^{(i)})\\
&=\sum\limits_{i=1}\limits^N\log\pi_iP(X,z_1=i|\theta^{(i)})，
(有\sum\limits_{i=1}\limits^N\pi_i=1)\\
l(\pi,\eta)&=\sum\limits_{i=1}\limits^N\log\pi_iP(X,z_1=i|\theta^{(i)})+\eta(\sum\limits_{i=1}\limits^N\pi_i-1)\\
令\frac{\partial l}{\partial\pi_i}=0，&则P(X,z_1=i|\theta^{(i)})+\pi_i\eta=0\\
&\Rightarrow \sum\limits_{i=1}\limits^N(P(X,z_1=i|\theta^{(i)})+\pi_i\eta)=0\\
&\Rightarrow \eta=-P(X|\theta^{(i)})，\pi_i=\frac{P(X,z_1=i|\theta^{(i)})}{P(X|\theta^{(i)})}
\end{align*}
$$

---

$$
\begin{align*}
\sum\limits_Z P(X,Z|\theta^{(i)})\sum\limits_{t=1}\limits^{T-1}\log a_{z_tz_{t+1}}
&=\sum\limits_{z_1}\dots\sum\limits_{z_T}P(X,z_1,\dots,z_T|\theta^{(i)})\sum\limits_{t=1}\limits^{T-1}\log a_{z_tz_{t+1}}\\
&=\sum\limits_{t=1}\limits^{T-1}\log a_{z_tz_{t+1}}\sum\limits_{i=1}\limits^N\sum\limits_{j=1}\limits^NP(X,z_t=i,z_{t+1}=j|\theta^{(i)})\\
&=\sum\limits_{i=1}\limits^{N}\sum\limits_{j=1}\limits^{N}\sum\limits_{t=1}\limits^{T-1}\log a_{ij}P(X,z_t=i,z_{t+1}=j|\theta^{(i)})，(有\sum\limits_{j=1}\limits^Na_{ij}=1)\\
l(a_{ij},\eta)&=\sum\limits_{i=1}\limits^{N}\sum\limits_{j=1}\limits^{N}\sum\limits_{t=1}\limits^{T-1}\log a_{ij}P(X,z_t=i,z_{t+1}=j|\theta^{(i)})+\eta(\sum\limits_{i=1}\limits^Na_{ij}-1)\\
令\frac{\partial l}{\partial a_{ij}}=0,&则\sum\limits_{t=1}\limits^{T-1}P(X,z_t=i,z_{t+1}=j|\theta^{(i)})+a_{ij}\eta=0\\
&\Rightarrow \sum\limits_{j=1}\limits^{N}(\sum\limits_{t=1}\limits^{T-1}P(X,z_t=i,z_{t+1}=j|\theta^{(i)})+a_{ij}\eta)=0\\
&\Rightarrow\eta=-\sum\limits_{t=1}\limits^{T-1}P(X,z_t=i|\theta^{(i)})，a_{ij}=\frac{\sum\limits_{t=1}\limits^{T-1}P(X,z_t=i,z_{t+1}=j|\theta^{(i)})}{\sum\limits_{t=1}\limits^{T-1}P(X,z_t=i|\theta^{(i)})}
\end{align*}
$$

---

$$
\begin{align*}
\sum\limits_Z P(X,Z|\theta^{(i)})\sum\limits_{t=1}\limits^{T}\log b_{z_t}(x_t)
&=\sum\limits_{z_1}\dots\sum\limits_{z_T}P(X,z_1,\dots,z_T|\theta^{(i)})\sum\limits_{t=1}\limits^{T}\log b_{z_t}(x_t)\\
&=\sum\limits_{t=1}\limits^{T}\log b_{z_t}(x_t)\sum\limits_{j=1}\limits^NP(X,z_t=j|\theta^{(i)})\\
&=\sum\limits_{j=1}\limits^N\sum\limits_{t=1}\limits^{T}\log b_{z_t}(x_t)P(X,z_t=j|\theta^{(i)})，(有\sum\limits_{k=1}\limits^Mb_j(k)=1)\\
l(b_j(k),\eta)&=\sum\limits_{j=1}\limits^N\sum\limits_{t=1}\limits^{T}\log b_{z_t}(x_t)P(X,z_t=j|\theta^{(i)})+\eta(\sum\limits_{k-1}\limits^Mb_j(k)-1)\\
令\frac{\partial l}{\partial b_j(k)}=0，&则\sum\limits_{t=1}\limits^{T}P(X,z_t=j,x_t=k|\theta^{(i)})+b_j(k)\eta=0\\
&\Rightarrow\sum\limits_{k=1}\limits^M(\sum\limits_{t=1}\limits^{T}P(X,z_t=j,x_t=k|\theta^{(i)})+b_j(k)\eta)=0\\
&\Rightarrow\eta=-\sum\limits_{t=1}\limits^TP(X,z_t=j|\theta^{(i)})，b_j(k)=\frac{\sum\limits_{t=1}\limits^TP(X,z_t=j,x_t=k|\theta^{(i)})}{\sum\limits_{t=1}\limits^TP(X,z_t=j|\theta^{(i)})}
\end{align*}
$$

**(d)**
$$
\begin{align*}
\gamma_t(i)&=P(z_t=q_i|X,\theta)
=\frac{P(z_t=q_i,X|\theta)}{P(X|\theta)}
=\frac{\alpha_t(i)\beta_t(i)}{\sum\limits_{j=1}\limits^N\alpha_t(j)\beta_t(j)}\\
\xi_t(i,j)&=P(z_t=q_i,z_{t+1}=q_j|X,\theta)
=\frac{P(z_t=q_i,z_{t+1}=q_j,X|\theta)}{P(X|\theta)}\\
&=\frac{\alpha_t(i)a_{ij}b_j(x_{t+1})\beta_{t+1}(j)}{\sum\limits_{i=1}\limits^N\sum\limits_{j=1}\limits^N\alpha_t(i)a_{ij}b_j(x_{t+1})\beta_{t+1}(j)}
\end{align*}
$$
**(e)**
$$
\begin{align*}
\pi_i&=\gamma_1(i)，
a_{ij}=\frac{\sum\limits_{t=1}\limits^{T-1}\xi_t(i,j)}{\sum\limits_{t=1}\limits^{T-1}\gamma_t(i)}，
b_j(k)=\frac{\sum\limits_{t=1,x_i=k}\limits^{T-1}\gamma_t(j)}{\sum\limits_{t=1}\limits^T\gamma_t(j)}
\end{align*}
$$
**伪代码：**

```pseudocode
init theta = (pi(0), A(0), B(0)), n = 0
while n < max and theta is not best
{
	a_ij(n+1) = sum(xi_t[i][j])/sum(gamma_t[i])
	b_j[k](n+1) = sum(gamma_t[j])/sum(gamma_t[j])
	pi_i(n+1) = gamma_1(i)
	n = n + 1
}
return theta = (pi, A, B)
```



## 3 LDA 实现

