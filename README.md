# OptStopRandNNCpp
 a simple C++ code for "Optimal Stopping via Randomized Neural Networks"
 
 Original Paper: https://arxiv.org/abs/2104.13669
 
 Their official code: https://github.com/HeKrRuTe/OptStopRandNN
 
 

used library: 
- Eigen(https://eigen.tuxfamily.org/index.php?title=Main_Page): a C++ template library for linear algebra
- EigenRand(https://github.com/bab2min/EigenRand): a random distribution generator for Eigen

구현하기에 매우 편리할 정도로 알고리즘이 단순해서... 구현해 보았다.

대략적으로 논문에서 제시한 것과 비슷한 결과가 나오는데,
- RFQI는 정직하게 inverse matrix를 계산하여 곱해 주기 때문인지 계산 시간이 오래 걸린다. RLSM, LSM과 같이 svd를 이용한 solver 를 사용하면 개선이 될 것 같다.
- RLSM은 C++을 사용해서인지 논문보다 계산 속도가 매우 빠르다.
- LSM은 차원이 낮을 때는 RLSM보다도 빠르다. dimension 3이하일 때는 RLSM조차도 쓸 이유가 없어 보인다. 그러나, 결과값은 그래도 RLSM이 더 높고 (물론 RFQI가 제일 높다), 차원이 증가하면 LSM의 소요시간은 매우 빠르게 증가한다.

LSM은 basis function을 그냥 일반 2차 다항식으로 해서인지 논문보다 결과가 항상 나쁘다. 논문에서는 Laguerre polynomial을 사용하였다. 

사용하기에 따라서는 RLSM, RFQI 모두 충분히 실제로 사용할 가치는 있어 보인다.
