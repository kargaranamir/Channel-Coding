clear

G=[1,0,0,1,1,1,0;0,1,0,0,1,1,1;0,0,1,1,1,0,1];

H=gen2par(G);
syndt=[
       1 1 1 0 1 0 0 0 0 0 0;
       0 1 1 1 0 1 0 0 0 0 0;
       1 1 0 1 0 0 1 0 0 0 0;
       1 0 0 0 0 0 0 1 0 0 0;
       0 1 0 0 0 0 0 0 1 0 0;
       0 0 1 0 0 0 0 0 0 1 0;
       0 0 0 1 0 0 0 0 0 0 1;
       0 0 1 1 0 0 0 0 0 1 1;
       0 1 1 0 0 0 0 0 1 1 0;
       1 1 0 0 0 0 0 1 1 0 0;
       0 1 0 1 0 0 1 1 0 0 0;
       1 0 1 0 0 1 1 0 0 0 0;
       1 0 0 1 1 1 0 0 1 0 0;
       1 1 1 1 0 1 0 1 0 0 0;
       1 0 1 1 0 1 1 0 0 0 1;
       0 0 0 0 0 0 0 0 0 0 0
        ];
v=[ 
    0 0 0  0 0 0 0 0 0 0;
    0 0 1  0 0 1 1 1 0 1;
    0 1 0  0 1 0 0 1 1 1;
    0 1 1  0 1 1 1 0 1 0;
    1 0 0  1 0 0 1 1 1 0;
    1 0 1  1 0 1 0 0 1 1;
    1 1 0  1 1 0 1 0 0 1;
    1 1 1  1 1 1 0 1 0 0;
   ];

data_uncoded= randi([0 1],1,10000);

data_uncoded_reshape=reshape([data_uncoded 0 0],[],3);

codewords=mod(data_uncoded_reshape*G,2);
codewords_reshape_padded=[reshape(codewords,1,[]) zeros(1,2)];

data_coded_2=bin2dec(int2str(reshape(codewords_reshape_padded,[],1)));
data_coded_8=bin2dec(int2str(reshape(codewords_reshape_padded,[],3)));
data_coded_16=bin2dec(int2str(reshape(codewords_reshape_padded,[],4)));

data_coded_bpsk=pskmod(data_coded_2,2);
data_coded_8psk=pskmod(data_coded_8,8);
data_coded_16psk=pskmod(data_coded_16,16);

data_uncoded_8=bin2dec(int2str(reshape([data_uncoded zeros(1,8)],[],3)));
data_uncoded_16=bin2dec(int2str(reshape([data_uncoded zeros(1,8)],[],4)));

data_uncoded_bpsk=pskmod(data_uncoded,2);
data_uncoded_8psk=pskmod(data_uncoded_8,8);
data_uncoded_16psk=pskmod(data_uncoded_16,16);

rub=[];
ru8=[];
ru16=[];
rcb=[];
rc8=[];
rc16=[];
parfor snr_db=-10:1:20
    
    data_coded_bpsk_AWGN_dmod=reshape(decimalToBinaryVector(pskdemod(awgn(data_coded_bpsk,snr_db*7/4),2)),1,[]);
    data_coded_8psk_AWGN_dmod=reshape(decimalToBinaryVector(pskdemod(awgn(data_coded_8psk,snr_db*7/4),8),3),1,[]);
    data_coded_16psk_AWGN_dmod=reshape(decimalToBinaryVector(pskdemod(awgn(data_coded_16psk,snr_db*7/4),16),4),1,[]);
    
    data_uncoded_bpsk_AWGN_dmod=reshape(decimalToBinaryVector(pskdemod(awgn(data_uncoded_bpsk',snr_db),2)),1,[]);
    data_uncoded_8psk_AWGN_dmod=reshape(decimalToBinaryVector(pskdemod(awgn(data_uncoded_8psk,snr_db),8),3),1,[]);
    data_uncoded_16psk_AWGN_dmod=reshape(decimalToBinaryVector(pskdemod(awgn(data_uncoded_16psk,snr_db),16),4),1,[]);
    
    data_coded_bpsk_reshape=reshape(data_coded_bpsk_AWGN_dmod(1:length(data_coded_bpsk_AWGN_dmod)-2),7,[]);
    data_coded_8psk_reshape=reshape(data_coded_8psk_AWGN_dmod(1:length(data_coded_8psk_AWGN_dmod)-2),7,[]);
    data_coded_16psk_reshape=reshape(data_coded_16psk_AWGN_dmod(1:length(data_coded_16psk_AWGN_dmod)-2),7,[]);
    
    data_coded_bpsk_synd=mod(H*data_coded_bpsk_reshape,2);
    data_coded_8psk_synd=mod(H*data_coded_8psk_reshape,2);
    data_coded_16psk_synd=mod(H*data_coded_16psk_reshape,2);
    
    epaternb=[];
    epatern8=[];
    epatern16=[];
    
   for i=1:length(data_coded_bpsk_synd(1,:))
       temp1=data_coded_bpsk_synd(:,i);
       temp2=data_coded_8psk_synd(:,i);
       temp3=data_coded_16psk_synd(:,i);
       for j=1:length(syndt(:,1))
           if (temp1==syndt(j,1:4)')
                epaternb=[epaternb syndt(j,5:end)'];
           end
           if (temp2==syndt(j,1:4)')
                epatern8=[epatern8 syndt(j,5:end)'];
           end
           if (temp3==syndt(j,1:4)')
                epatern16=[epatern16 syndt(j,5:end)'];
           end
           
       end
   end
   
   data_coded_bpsk_noerr=mod(data_coded_bpsk_reshape+epaternb,2);
   data_coded_8psk_noerr=mod(data_coded_8psk_reshape+epatern8,2);
   data_coded_16psk_noerr=mod(data_coded_16psk_reshape+epatern16,2);
   data_coded_bpsk_dcdr=[];
   data_coded_8psk_dcdr=[];
   data_coded_16psk_dcdr=[];
    for i=1:length(data_coded_bpsk_noerr(1,:))
       temp1=data_coded_bpsk_noerr(:,i);
       temp2=data_coded_8psk_noerr(:,i);
       temp3=data_coded_16psk_noerr(:,i);
       c1=0;
       c2=0;
       c3=0;
       for j=1:length(v(:,1))
           if (temp1==v(j,4:end)')
                data_coded_bpsk_dcdr=[data_coded_bpsk_dcdr v(j,1:3)'];
           else 
                c1=c1+1;
                if (c1==(length(v(:,1))))
                    data_coded_bpsk_dcdr=[data_coded_bpsk_dcdr 2*ones(3,1)];
                end
           end
           if (temp2==v(j,4:end)')
                data_coded_8psk_dcdr=[data_coded_8psk_dcdr v(j,1:3)'];
           else  
                c2=c2+1;
                if (c2==(length(v(:,1))))
                    data_coded_8psk_dcdr=[data_coded_8psk_dcdr 2*ones(3,1)];
                end
           end
           if (temp3==v(j,4:end)')
                data_coded_16psk_dcdr=[data_coded_16psk_dcdr v(j,1:3)'];
           else   
                c3=c3+1;
                if (c3==(length(v(:,1))))
                    data_coded_16psk_dcdr=[data_coded_16psk_dcdr 2*ones(3,1)];
                end
           end
       end
   end  
    data_coded_bpsk_dcdr_rshape=reshape(data_coded_bpsk_dcdr,1,[]);
    data_coded_8psk_dcdr_rshape=reshape(data_coded_8psk_dcdr,1,[]);
    data_coded_16psk_dcdr_rshape=reshape(data_coded_16psk_dcdr,1,[]);
    
    data_coded_bpsk_dcdr_rshape=data_coded_bpsk_dcdr_rshape(1:(length(data_coded_bpsk_dcdr_rshape)-2));
    data_coded_8psk_dcdr_rshape=data_coded_8psk_dcdr_rshape(1:(length(data_coded_8psk_dcdr_rshape)-2));
    data_coded_16psk_dcdr_rshape=data_coded_16psk_dcdr_rshape(1:(length(data_coded_16psk_dcdr_rshape)-2));
    
    data_uncoded_bpsk_rec=data_uncoded_bpsk_AWGN_dmod;
    data_uncoded_8psk_rec=data_uncoded_8psk_AWGN_dmod(1:length(data_uncoded_8psk_AWGN_dmod)-8);
    data_uncoded_16psk_rec=data_uncoded_16psk_AWGN_dmod(1:length(data_uncoded_16psk_AWGN_dmod)-8);
    
    [ntemp,rtemp]=biterr(data_coded_bpsk_dcdr_rshape,data_uncoded);
    rcb=[rcb rtemp];
    [ntemp,rtemp]=biterr(data_coded_8psk_dcdr_rshape,data_uncoded);
    rc8=[rc8 rtemp];
    [ntemp,rtemp]=biterr(data_coded_16psk_dcdr_rshape,data_uncoded);
    rc16=[rc16 rtemp];
    [ntemp,rtemp]=biterr(data_uncoded_bpsk_rec,data_uncoded);
    rub=[rub rtemp];
    [ntemp,rtemp]=biterr(data_uncoded_8psk_rec,data_uncoded);
    ru8=[ru8 rtemp];
    [ntemp,rtemp]=biterr(data_uncoded_16psk_rec,data_uncoded);
    ru16=[ru16 rtemp];
    
end
snr_db=-10:1:20;
figure
semilogy(snr_db,rub,snr_db,rcb)
grid on
title('BPSK')
legend('Uncoded','Coded')
ylabel('bit error rate')
xlabel('SNR(dB)')
figure
semilogy(snr_db,ru8,snr_db,rc8)
grid on
title('8PSK')
legend('Uncoded','Coded')
ylabel('bit error rate')
xlabel('SNR(dB)')
figure
semilogy(snr_db,ru16,snr_db,rc16)
grid on
title('16PSK')
legend('Uncoded','Coded')
ylabel('bit error rate')
xlabel('SNR(dB)')