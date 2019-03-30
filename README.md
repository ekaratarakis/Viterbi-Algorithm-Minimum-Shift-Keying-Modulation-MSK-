# Viterbi-Algorithm-Minimum-Shift-Keying-Modulation-MSK-
This project is the implementation of the Viterbi algorithm for estimating the maximum aposteriori probability on the Minimum Shift Keying modulation (MSK). To begin with, it is important to study the Bit Error Rate (BER) of the MSK and to do so we had to implement the MSK modulation in Matlab by pairing the symbols and creating QPSK symbols from these pairs. The conversion from MSK to QPSK was based on the recursive equations below

![image convert from MSK to QPSK](https://github.com/ekaratarakis/Viterbi-Algorithm-Minimum-Shift-Keying-Modulation-MSK-/blob/master/Images/convert.jpg)

Having defined the Signal to Noise Ration (SNR) like in the equation below 

![image SNR](https://github.com/ekaratarakis/Viterbi-Algorithm-Minimum-Shift-Keying-Modulation-MSK-/blob/master/Images/snr.jpg)

and by giving **SNR = 5** was simulated the equation

![image equation 11](https://github.com/ekaratarakis/Viterbi-Algorithm-Minimum-Shift-Keying-Modulation-MSK-/blob/master/Images/eq11.jpg) 

and estimated the BER as **BER = 0.075249**

Subsequently, we repeat the above procedure for SNR = 6,7,8,9,10,11,12  
## Viterbi Algorithm
Now, lets go deeper into the juicy part and implement the Viterbi algorithm. The equivalent baseband signal can be expressed with the equations below

![image equivalent baseband signal](https://github.com/ekaratarakis/Viterbi-Algorithm-Minimum-Shift-Keying-Modulation-MSK-/blob/master/Images/viterbi.jpg)

which allows us to simulate the equation 

![image equation 11](https://github.com/ekaratarakis/Viterbi-Algorithm-Minimum-Shift-Keying-Modulation-MSK-/blob/master/Images/eq11.jpg)

for a quite large number of symbols and implement the Viterbi algorithm to obtain the maximum likelihood sequence of the symbols. In order to implement that we had to define the constant vectors s^{-1} and s^{1} which refer to receive symbol xn = -1 or xn = +1 respectively. In order to calculate the weights in every transition we used the equation 

![image weights](https://github.com/ekaratarakis/Viterbi-Algorithm-Minimum-Shift-Keying-Modulation-MSK-/blob/master/Images/weights.jpg)

Additionaly, in order to find the maximum likelihood sequence we needed to use the trellis diagram below

![image trellis](https://github.com/ekaratarakis/Viterbi-Algorithm-Minimum-Shift-Keying-Modulation-MSK-/blob/master/Images/trellis.jpg)
