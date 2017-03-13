pro movieindex,fileout,files, titlehtml=titlehtml

if (n_elements(titlehtml) eq 0) then titlehtml='Javascript movie'

movie0 = $
'<HTML>'+string('0A'XB)  $
+''+string('0A'XB) $
+'<HEAD>'+string('0A'XB) $
+'<TITLE>' + titlehtml +'</TITLE>'+string('0A'XB) $  ;;;;; <--- title
+'</HEAD>'+string('0A'XB) $
+''+string('0A'XB) $
+'<BODY BGCOLOR="White">'+string('0A'XB) $
+''+string('0A'XB) $
+'<CENTER>'+string('0A'XB) $
+''+string('0A'XB) $
+'<H1> '+titlehtml+' </H1>'+string('0A'XB) $ ;;;;; <--- title
+'<P>'+string('0A'XB) $
+''+string('0A'XB) $
+'<TABLE BORDER="10" CELLPADDING="8">'+string('0A'XB) $
+'<TR>'+string('0A'XB) $
+'<TD align="center"> '+string('0A'XB) 

movie1 = $
'</TD></TR>'+string('0A'XB) $
+'</TABLE>'+string('0A'XB) $
+'<P>'+string('0A'XB) $
+''+string('0A'XB) $
+'<FORM NAME="form">'+string('0A'XB) $
+' <INPUT TYPE=button VALUE="Start" onClick="start_play();">'+string('0A'XB) $
+' <INPUT TYPE=button VALUE="Pause" onClick="stop_play();">'+string('0A'XB) $
+' <INPUT TYPE=button VALUE="Faster" onClick="delay/=inc; show_delay();">'+string('0A'XB) $
+' <INPUT TYPE=button VALUE="Slower" onClick="delay*=inc; show_delay();">'+string('0A'XB) $
+' <INPUT TYPE=button VALUE="Step" onClick="step();">'+string('0A'XB) $
+' <INPUT TYPE=button NAME="direction" VALUE="Reverse" onClick="reverse();">'+string('0A'XB) $
+' <INPUT TYPE=button VALUE="Swing_Mode" onClick="swing_mode();">'+string('0A'XB) $
+' <BR>'+string('0A'XB) $
+' Frame: <INPUT TYPE=text VALUE="" NAME="frame" SIZE=30> '+string('0A'XB) $
+'<!-- &nbsp;Delay (ms): <INPUT TYPE=text VALUE="" NAME="delay" SIZE=6> //-->'+string('0A'XB) $
+''+string('0A'XB) $
+'</FORM>'+string('0A'XB) $
+'</CENTER>'+string('0A'XB) $
+''+string('0A'XB) $
+'<P>'+string('0A'XB) $
+'<HR>'+string('0A'XB) $
+'<B>Document</B>: <I><SCRIPT>document.write(document.title);</SCRIPT></I><BR>'+string('0A'XB) $
+'<B>URL</B>: <I><SCRIPT>document.write(document.URL);</SCRIPT></I><BR>'+string('0A'XB) $
+'<B>Last Update</B>: <I><SCRIPT>document.write(document.lastModified);</SCRIPT></I><BR>'+string('0A'XB) $
+''+string('0A'XB) $
+'</BODY>'+string('0A'XB) $
+'</HTML>'+string('0A'XB) $
+''+string('0A'XB) $
+'<SCRIPT>'+string('0A'XB) $
+'<!--'+string('0A'XB) $
+''+string('0A'XB) $
+'// Javascript program to produce fast animations by reading from cache'+string('0A'XB) $
+'// Written: 7 July 1997, Zarro, NASA/GSFC'+string('0A'XB) 

movie2 = $
'var num_loaded_images = 0;'+string('0A'XB) $
+'var frame=-1;'+string('0A'XB) $
+'var timeout_id=null;'+string('0A'XB) $
+'var dir=1, playing=0;'+string('0A'XB) $
+'var run = 0;'+string('0A'XB) $
+'var bname = "Reverse";'+string('0A'XB) $
+'var swingon = 0;'+string('0A'XB) $
+''+string('0A'XB) $
+'// function to count images as they are loaded into cache'+string('0A'XB) $
+''+string('0A'XB) $
+'function count_images() '+string('0A'XB) $
+'{ '+string('0A'XB) $
+' if (++num_loaded_images == imax) '+string('0A'XB) $
+'  animate(); '+string('0A'XB) $
+' else {'+string('0A'XB) $
+'  document.animation.src=images[num_loaded_images-1].src;'+string('0A'XB) $
+'  document.form.frame.value="Loading "+num_loaded_images+" of "+imax; '+string('0A'XB) $
+' }'+string('0A'XB) $
+'}'+string('0A'XB) $
+''+string('0A'XB) $
+'temp = new Array(imax);          '+string('0A'XB) 


movie3 = $
+''+string('0A'XB) $
+'// actual loading is done here'+string('0A'XB) $
+' '+string('0A'XB) $
+'images = new Array(imax);'+string('0A'XB) $
+'for (var i = 0 ; i < imax; i++) {'+string('0A'XB) $
+' images[i]= new Image();'+string('0A'XB) $
+' images[i].onload=count_images;'+string('0A'XB) $
+' images[i].src=temp[i];'+string('0A'XB) $
+'}'+string('0A'XB) $
+''+string('0A'XB) $
+'// function to start movie'+string('0A'XB) $
+''+string('0A'XB) $
+'function start_play() {'+string('0A'XB) $
+' if (playing == 0) {'+string('0A'XB) $
+'  if (timeout_id == null && num_loaded_images==imax) animate();'+string('0A'XB) $
+' }'+string('0A'XB) $
+'} '+string('0A'XB) $
+''+string('0A'XB) $
+'// function to stop movie'+string('0A'XB) $
+''+string('0A'XB) $
+'function stop_play() {'+string('0A'XB) $
+' if (timeout_id) clearTimeout(timeout_id); '+string('0A'XB) $
+'  timeout_id=null;'+string('0A'XB) $
+'  playing = 0;'+string('0A'XB) $
+'} '+string('0A'XB) $
+''+string('0A'XB) $
+'// function to control swing mode'+string('0A'XB) $
+'function swing_mode() {'+string('0A'XB) $
+'  if (swingon) swingon=0; else swingon=1;'+string('0A'XB) $
+'}'+string('0A'XB) $
+'  '+string('0A'XB) $
+''+string('0A'XB) $
+'// function to do the animation when all images are loaded'+string('0A'XB) $
+''+string('0A'XB) $
+'function animate()'+string('0A'XB) $
+'{'+string('0A'XB) $
+' var j;'+string('0A'XB) $
+' frame=(frame+dir+imax)%imax;'+string('0A'XB) $
+' j=frame+1;'+string('0A'XB) $
+' document.animation.src=images[frame].src;'+string('0A'XB) $
+' document.form.frame.value="Displaying "+j+" of "+imax;'+string('0A'XB) $
+' if (swingon && (j == imax || frame == 0)) reverse();'+string('0A'XB) $
+' timeout_id=setTimeout("animate()",delay);'+string('0A'XB) $
+' playing=1;'+string('0A'XB) $
+'}'+string('0A'XB) $
+''+string('0A'XB) $
+''+string('0A'XB) $
+'// function to control stepping thru each frame'+string('0A'XB) $
+''+string('0A'XB) $
+'function step()'+string('0A'XB) $
+'{'+string('0A'XB) $
+' var j;'+string('0A'XB) $
+' if (timeout_id) clearTimeout(timeout_id); timeout_id=null;'+string('0A'XB) $
+' frame=(frame+dir+imax)%imax;'+string('0A'XB) $
+' j=frame+1;'+string('0A'XB) $
+' document.animation.src=images[frame].src;'+string('0A'XB) $
+' document.form.frame.value="Displaying "+j+" of "+imax;'+string('0A'XB) $
+' playing=0;'+string('0A'XB) $
+'}'+string('0A'XB) $
+''+string('0A'XB) $
+'// function to control direction of animation'+string('0A'XB) $
+''+string('0A'XB) $
+'function reverse()'+string('0A'XB) $
+'{'+string('0A'XB) $
+' dir=-dir;'+string('0A'XB) $
+' if (dir > 0) document.form.direction.value="Reverse"; bname="Reverse";'+string('0A'XB) $
+' if (dir < 0) document.form.direction.value="Forward"; bname="Forward";'+string('0A'XB) $
+'}'+string('0A'XB) $
+''+string('0A'XB) $
+'// function to display delay between frames (not implemented yet)'+string('0A'XB) $
+''+string('0A'XB) $
+'function show_delay()'+string('0A'XB) $
+'{'+string('0A'XB) $
+'// document.form.delay.value=delay;'+string('0A'XB) $
+'}'+string('0A'XB) $
+''+string('0A'XB) $
+'// -->'+string('0A'XB) $
+'</SCRIPT>'+string('0A'XB)

mx=n_elements(files)

openw,lun,/get_lun,fileout

printf,lun,movie0
printf,lun,'<img src = "'+files[0]+'" NAME=animation>'+string('0A'XB)
printf,lun,movie1
printf,lun,'var imax = '+strtrim(mx,2)+', inc = 1.50, delay = 250;'
printf,lun,movie2
for m=0,mx-1 do begin
printf,lun,'temp['+strtrim(m,2)+']="'+files[m]+'";'
endfor
printf,lun,movie3

close,lun
free_lun,lun
end
