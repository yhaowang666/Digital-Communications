clc;clear all;close all;
N = 1000000;
s = source(N); %��Դ���������и���ΪN

b = hamming_encoding(s);
y = hamming_decoding(b);