       program junk

       double precision a


       a = float(30)/float(140)
       print *,sizeof(float(30)),a
       a = real(30)/real(140)
       print *,sizeof(real(30)),a
       a = dble(30)/dble(140)
       print *,sizeof(dble(30)),a
       stop
       end
