# Experiment problems of SDC

## 2021.11.2

### 1. Overview

My implementation involves more error and orthogonality loss compare with MATLAB function eig() when the block size is big, and when I set block size as 1, the boost of quality of the solution does improve but just doing slightly better than eig(). In contrast, Yuji's official code act just like what he said in his paper, even I set block size 32 or more. I try to improve my implementation latter on. The followings are some of the improvements have been done.

#### 1) Orthogonality of QDWH

I found that the unitary matrix computed by my QDWH has worse orthogonality than Yuji's,  I found the terminate condition of Yuji's code involving an equality
$$
|1 - L| < \epsilon
$$
<span style = "color:red">which you need to check the related paper to figure it out later.</span> The worst thing is, even I apply the above change to my code, the orthogonality of my QDWH cannot reach Yuji's (but it has considerably improvement), <span style = "color:red"> it still needs time to figure out why. </span>

#### 2) Subspace iteration

After some improvement of QDWH, the quality of eigenvalue/eigenvector of my QDWH still struggle, <span style = "color:red">so next time I think I need to carefully check this part for better performance.</span>

## 2021.11.17

I finally found the above problems are caused by Newton-Schulz postprocessing, it goes as follows
$$
U = \frac{3}{2} U - \frac{1}{2}U(U'U)
$$
I gonna searching for this...

So the next target is to find out QDWH and NS postprocessing!

### 2021.12.7

Some notes for calling LAPACK...

You should call DGETRF before DGECON to finish LU factorization for work in DGECON...

I'm stupid lol... Trying to figure out why my DGECON computes such a big estimation for 1-Norm condition number...

### 2021.12.10

Finished v1.1 branch, try symmetric operation inside the QDWH algorithm. But the performance doesn't look really good... It is not able to beat dsyev() in small cases, and in large cases, it is faster than dsyev() but, I haven't implement other part of the QDWH-eig. I tried some optimization but it doesn't work very well. The biggest problem might be the memory copy.

But still great work today, keep calm and carry on!

### 2021.12.11

Try more symmetric operation and happy to see that the performance improves! Next goal is to implement the cholesky version of QDWH for well conditioned U.

### 2021.12.14

Right now I have finished the Cholesky part and see a great boost of the performance, keep on going!
Also, I delete Uprev and store the polar factor in the input matrix A (upper triangular)
