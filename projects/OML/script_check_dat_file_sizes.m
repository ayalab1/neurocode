q = dir('amplifier.dat');
sizes(1) = q.bytes;
q = dir('digitalin.dat');
sizes(2) = q.bytes;
q = dir('analogin.dat');
sizes(3) = q.bytes;
q = dir('time.dat');
sizes(4) = q.bytes;
q = dir('auxiliary.dat');
sizes(5) = q.bytes;

% sizes./sizes(3)
any(rem(sizes/(sizes(2)),1))