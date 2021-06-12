# OptStopRandNNCpp
 a simple C++ code for "Optimal Stopping via Randomized Neural Networks"
 
 Original Paper: https://arxiv.org/abs/2104.13669
 
 Their official code: https://github.com/HeKrRuTe/OptStopRandNN
 
 

used library: 
- Eigen(https://eigen.tuxfamily.org/index.php?title=Main_Page): a C++ template library for linear algebra
- EigenRand(https://github.com/bab2min/EigenRand): a random distribution generator for Eigen

일단은 non-Markovian case에는 별로 관심이 없어서 RLSM, RFQI만 C++로 구현을 해 보았다.
RFQI는 그럭저럭 결과가 잘 나오는 듯하지만 RLSM은 논문에서 말하는 정도의 결과가 나오지 않는다.
inverse 를 구하는 행렬이 ill-posed 라서 그런 것 같아서, Eigen에서 제공하는 SVD 함수를 이용해 least square solution을 구해 보았지만, 여전히 안정적이지는 않은 모양이다.

python 코드로는 이 문제를 어떻게 해결했는지 궁금하구먼...
