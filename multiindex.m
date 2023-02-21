function [ ind ] = multiindex( d,k )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if d==1 
    ind(1)=k;
else
    nbot=1;
    for i=1:k-d+1
        indextemp=multiindex(d-1,k-i);
   
     [s,t]=size(indextemp);
        ntop=nbot+s-1;
        ind(nbot:ntop,1)=i*ones(s,1);
        ind(nbot:ntop,2:d)=indextemp;
        nbot=ntop+1;
    end
end

end
