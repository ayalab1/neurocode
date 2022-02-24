count=0;
mode=4;
while true 

                    % ask arduino send info
                    flushinput(aHandle);
                    flushoutput(aHandle);
                    fwrite(aHandle,4,'uchar');
                    if(aHandle.BytesAvailable<6) % receive it if there is one
                        ontinue
                    end
                    pause(0.1);
                    b = fread(aHandle,6);
                    ArduinoCheckPort = typecast(uint8(b(1)),'int8'); 
                    plot(count,b);
                    pause(1);
                    count=count+1 ;
                    mode=4;
                    
                   
                  
end
                  