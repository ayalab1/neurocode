function out = sec2min(sec)

out = [floor(sec/60) rem(sec,60)];
out = [num2str(out(1)) 'm ' num2str(round(out(2))) 's'];