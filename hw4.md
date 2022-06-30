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
<div style="page-break-after: always;"></div>

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

<div style="page-break-after: always;"></div>

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
	a_ij(n+1) = sum(xi_t[i][j])/sum(gamma_t[i])
	b_j[k](n+1) = sum(gamma_t[j])/sum(gamma_t[j])
	pi_i(n+1) = gamma_1(i)
	n = n + 1
return theta = (pi, A, B)
```



## 3 LDA 实现

**(a)**

```pseudocode
init params
while epoch < max_epochs and delta > tolerance
	while d < D
		while i < max_iter_d and delta >= tol_d
			for n = 0 to n_word_types
				compute phi
		compute gamma
		update bound
	update log_betas
	update alpha
	compute delta
	store new value of bound
	output
```

**(b)** 见 main.py

**(c)**

K = 5

```
0 孩子 学生 医院 学校 老师 家长 医生 告诉
1 司机 车辆 女士 交警 发生 乘客 现场 事故
2 民警 男子 警方 嫌疑人 发现 派出所 犯罪 报警
3 公司 万元 发现 银行 工作人员 快递 相关 工作
4 法院 老人 儿子 母亲 发现 父亲 女儿 死亡
```

K = 10

```
0 学生 学校 老师 网友 同学 家长 女生 微博
1 医院 医生 治疗 手术 患者 检查 情况 病人
2 民警 警方 嫌疑人 男子 手机 犯罪 发现 李某
3 男子 司机 民警 车辆 发生 交警 现场 事故
4 警方 民警 男子 发现 女子 老人 派出所 调查
5 公司 万元 银行 发现 女士 工作人员 电话 快递
6 游客 李桂英 村民 发现 动物 保护 多年 村里
7 孩子 儿子 家长 妈妈 父母 女儿 父亲 家里
8 法院 赔偿 母亲 父亲 判决 女儿 离婚 妻子
9 工作 中国 时间 报道 美国 社会 公司 生活
```

K = 20

```
0 孩子 家长 女儿 幼儿园 女士 儿子 妈妈 老师
1 离婚 结婚 丈夫 两人 妻子 生活 男友 父母
2 儿子 父亲 母亲 孩子 父母 家里 女儿 家人
3 万元 公司 美国 银行 女性 工作 女士 介绍
4 公司 法院 赔偿 女士 万元 手机 被告 证明
5 医院 医生 患者 价格 发现 治疗 检查 市场
6 司机 车辆 交警 乘客 公交车 现场 发生 男子
7 嫌疑人 犯罪 警方 公安局 视频 案件 民警 调查
8 学生 学校 老师 家长 同学 网友 孩子 女生
9 游客 导游 旅游 旅行社 景区 刘先生 发现 标准
10 民警 男子 派出所 发现 警方 交警 处罚 旅客
11 法院 房屋 被告人 判处 判决 有期徒刑 证据 拆迁
12 小区 女子 业主 物业 居民 发现 保安 现场
13 医院 医生 手术 男子 治疗 家属 患者 发生
14 报道 动物 保护 英国 一只 中国 野生动物 发现
15 男子 民警 警方 李某 王某 妻子 发现 嫌疑人
16 警方 工作 发现 调查 人员 电话 死亡 公司
17 事故 发生 车辆 救援 驾驶 责任 保险公司 现场
18 老人 大爷 师傅 儿子 发现 老伴 下午 老太
19 手机 民警 发现 警方 盗窃 嫌疑人 超市 男子
```

**(d)** 分类效果最好的 K = 20，因为数据集本身文本量大，主题数量较多，如果 K 取值过小就会因 topic 过少导致分类效果不佳，选用较大的 K 分类更为精准细致。
