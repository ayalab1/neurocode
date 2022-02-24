function varargout = CountBar(counts,portion)

% Each row sums up to a total. Row elements (diff columns) compared to each other
hold off;
sums = counts(:,end);
counts(:,end) = [];
[nRows nCols] = size(counts);

if exist('portion','var')>0 && portion==1
    ydata = counts./repmat(sums,1,nCols); ydata = ydata*100;
else
    ydata = counts;
end
ydata = ydata;
b = bar(ydata); hold on;

vers = version;
if num2str(vers(1))<8 % Code for Matlab 2010
	bars = get(b,'children'); %where 1st cell is first bars for each row
	for col = 1:nCols,
		xdata(:,col) = nanmean(get(bars{col},'xdata'));
	end
else
	for col = 1:nCols,
		xdata(:,col) = b(col).XData + b(col).XOffset;
	end
end
z = nan(nCols,nCols,nRows);

smallnumber = nanstd(ydata(:))/20;
smallnumberx = nanmean(diff(xdata(1,:)))/10; if isnan(smallnumberx), smallnumberx=0; end

for row = 1:nRows
    for i=1:nCols
        for j=1:row-1
            z(i,j,row) = zBinomialComparison(counts(row,i),sums(row),counts(j,i),sums(j));
            if abs(z(i,j,row))>1.96, %p<0.025 i.e. 5% two-tailed,
                % draw lines showing which columns are different:
                x1 = xdata(row,i) + abs(row-j)*smallnumberx;
                x2 = xdata(j,i) + abs(row-j)*smallnumberx;
                [~,this] = max(-counts([row j],i));
                if this==1, this=j; else, this=row; end;
%                 ydata(row,this) = ydata(row,this)+smallnumber;
                y1 = ydata(this,i) + smallnumber*abs(row-j);
                plot([x1 x2], [y1+smallnumber*sign(y1) y1+smallnumber*sign(y1)], 'k');
                plot([x1 x1], [y1+smallnumber*sign(y1) y1], 'k');
                plot([x2 x2], [y1+smallnumber*sign(y1) y1], 'k');
                % draw stars for the appropriate significance level:
                if abs(z(i,j,row))>3.3, %p<0.0005 i.e. 0.1% two-tailed
                    text(mean([x1 x2]), y1+smallnumber*sign(y1)*1.5, '***','HorizontalAlignment','center');
                elseif abs(z(i,j,row))>2.58, %p<0.005 i.e. 1% two-tailed
                    text(mean([x1 x2]), y1+smallnumber*sign(y1)*1.5, '**','HorizontalAlignment','center');
                else
                    text(mean([x1 x2]), y1+smallnumber*sign(y1)*1.5, '*','HorizontalAlignment','center');
                end
            end
        end
    end
end

if nargout>0,
    varargout{1} = b;
    varargout{2} = z;
end

%% Legacy (old code):
% function varargout = CountBar(counts,portion)
% 
% % Each row sums up to a total. Row elements (diff columns) compared to each other
% hold off;
% [nRows nCols] = size(counts);
% sums = sum(counts,2);
% 
% if exist('portion','var')>0 && portion==1,
%     ydata = counts./repmat(sum(counts,2),1,nCols);
% else
%     ydata = counts;
% end
% 
% b = bar(ydata); hold on;
% 
% vers = version;
% if num2str(vers(1))<8 % Code for Matlab 2010
% 	bars = get(b,'children'); %where 1st cell is first bars for each row
% 	for col = 1:nCols,
% 		xdata(:,col) = nanmean(get(bars{col},'xdata'));
% 	end
% else
% 	for col = 1:nCols,
% 		xdata(:,col) = b(col).XData + b(col).XOffset;
% 	end
% end
% z = nan(nCols,nCols,nRows);
% 
% smallnumber = nanstd(ydata(:))/20;
% smallnumberx = nanmean(diff(xdata(1,:)))/10;
% for row = 1:nRows,
%     for i=1:nCols,
%         for j=1:i-1,
%             z(i,j,row) = zBinomialComparison(counts(row,i),sums(row),counts(row,j),sums(row));
%             if abs(z(i,j,row))>1.96, %p<0.025 i.e. 5% two-tailed,
%                 % draw lines showing which columns are different:
%                 x1 = xdata(row,i) + abs(i-j)*smallnumberx;
%                 x2 = xdata(row,j) + abs(i-j)*smallnumberx;
%                 [~,this] = max(counts(row,[i j]));
%                 if this==1, this=i; else, this=j; end;
% %                 ydata(row,this) = ydata(row,this)+smallnumber;
%                 y1 = ydata(row,this) + smallnumber*abs(i-j);
%                 plot([x1 x2], [y1+smallnumber*sign(y1) y1+smallnumber*sign(y1)], 'k');
%                 plot([x1 x1], [y1+smallnumber*sign(y1) y1], 'k');
%                 plot([x2 x2], [y1+smallnumber*sign(y1) y1], 'k');
%                 % draw stars for the appropriate significance level:
%                 if abs(z(i,j,row))>3.3, %p<0.0005 i.e. 0.1% two-tailed
%                     text(mean([x1 x2]), y1+smallnumber*sign(y1)*1.5, '***','HorizontalAlignment','center');
%                 elseif abs(z(i,j,row))>2.58, %p<0.005 i.e. 1% two-tailed
%                     text(mean([x1 x2]), y1+smallnumber*sign(y1)*1.5, '**','HorizontalAlignment','center');
%                 else
%                     text(mean([x1 x2]), y1+smallnumber*sign(y1)*1.5, '*','HorizontalAlignment','center');
%                 end
%             end
%         end
%     end
% end
% 
% if nargout>0,
%     varargout{1} = z;
% end