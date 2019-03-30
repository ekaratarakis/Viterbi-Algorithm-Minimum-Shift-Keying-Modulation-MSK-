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
