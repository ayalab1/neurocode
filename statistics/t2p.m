function p = t2p(t,df)

%t2p - Transform a t-statistic into a p-value

p = 2*tcdf(-abs(t),df);