G=[1,0,0,1,1,1,0;0,1,0,0,1,1,1;0,0,1,1,1,0,1]
H=gen2par(G)

syndt=syndtable(H)

data_uncoded= randi([0 1],1,9999);
data=encode(data_uncoded,7,3,'linear',G);

bpskmod_data= pskmod(data',2);
bpskmod_udata=pskmod(data_uncoded',2);

brate=[];
burate=[];

for snr=-10:1:20
    bch_data=awgn(bpskmod_data,snr);
    bch_udata=awgn(bpskmod_udata,snr);
    
    bpskdmod_data=pskdemod(bch_data,2);
    bdata=decode(bpskdmod_data',7,3,'linear',G);
    [n,r]=biterr(data_uncoded,bdata);
    brate=[brate r];
    
    bpskdmod_udata=pskdemod(bch_udata,2);
    [n,r]=biterr(data_uncoded,bpskdmod_udata');
    burate=[burate r];
end

semilogy([-10:1:20],brate,[-10:1:20],burate)
legend('BPSK-Coded','BPSK-Uncoded')
grid on
title('BPSK-AWGN bit error')
ylabel('bit error rate')
xlabel('SNR(dB)')
