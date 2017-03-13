u0='<ul>'
u1='</ul>'
d0='<li>'
a0='<a href="'
a1='">'
a2='</a>'
i0='<img src="'
i1='">'

dirroot='/scr/s61/yokoyama/cans'
cd,dirroot

header= $
'<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">' $
+'<html lang="en-US">' $
+'<head>' $
+'<meta http-equiv="Content-Type" CONTENT="text/html; charset=ISO-8859-1">' $
+'<meta http-equiv="Content-Style-Type" CONTENT="text/css">' $
+'<title>' $
+'CANS' $
+'</title>' $
+'</head>' $
+'<body style="margin-left:30pt;margin-right:30pt;margin-top:30pt">'

footer= '</body>'+'</html>'

openw,unit,/get_lun,'index.html'
print,unit,header

dir0=file_search('cans?d',/test_directory,/mark_directory)
md0x=n_elements(dir0)


for md0=0,md0x-1 do begin

  printf,unit,file_basename(dir0[md0])
  printf,unit,u0

  dir1=file_search(dir0[md0]+'md*',/test_directory,/mark_directory)
  md1x=n_elements(dir1)

  for md1=0,md1x-1 do begin

    printf,unit,d0+a0+dir1[md1]+'index.html'+a1+file_basename(dir1[md1])+a2
    if (file_test(dir1[md1]+'README')) then begin
        printf,unit,' : <a href="'+dir1[md1]+'README'+a1+'README'+a2
    endif
    if (file_test(dir1[md1]+'README.txt')) then begin
        printf,unit,' : <a href="'+dir1[md1]+'README.txt'+a1+'README.txt'+a2
    endif
    if (file_test(dir1[md1]+'Readme.pdf')) then begin
        printf,unit,' : <a type="application/pdf" href="'+dir1[md1]+'Readme.pdf'+a1+'Readme.pdf'+a2
    endif

    print,dir1[md1]+'index.html'
    openw,unitmd,/get_lun,dir1[md1]+'index.html'
    printf,unitmd,header
    printf,unitmd,'<h2>'+file_basename(dir1[md1])+'</h2>'
    if (file_test(dir1[md1]+'Readme.pdf')) then begin
        printf,unitmd,'<a type="application/pdf" href="'+'Readme.pdf'+a1+'Readme.pdf'+a2+'<br>'
    endif
    if (file_test(dir1[md1]+'README')) then begin
        printf,unitmd,'<a href="'+'README'+a1+'README'+a2+'<br>'
    endif
    if (file_test(dir1[md1]+'README.txt')) then begin
        printf,unitmd,'<a href="'+'README.txt'+a1+'README.txt'+a2+'<br>'
    endif

    dir2=file_search(dir1[md1]+'*',/test_directory,/mark_directory)
    md2x=n_elements(dir2)

    for md2=0,md2x-1 do begin
      printf,unitmd,'<h3>'+file_basename(dir2[md2])+'</h3>'

      if (file_test(dir2[md2]+'moviedt/movie.html')) then begin
        printf,unitmd,'<a href="'+file_basename(dir2[md2])+'/moviedt/movie.html'+a1+'Movie'+a2+'<br>'
      endif

      pngfile=file_search(dir2[md2]+'*.png')
      md3x=n_elements(pngfile)

      for md3=0,md3x-1 do begin
        printf,unitmd,i0+file_basename(dir2[md2])+'/'+file_basename(pngfile[md3])+i1
      endfor

    endfor
    printf,unitmd,footer
    close,unitmd
    free_lun,unitmd

  endfor
  printf,unit,u1

endfor

printf,unit,footer
close,unit
free_lun,unit

end
