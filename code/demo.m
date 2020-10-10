clc;clear all;close all;
N = 1000000;
s = source(N); %信源产生，序列个数为N

b = hamming_encoding(s);
y = hamming_decoding(b);