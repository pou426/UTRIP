

CANS nc �⥸�塼��

       ver.0  2001.9.7

----------------------------------------------------------------------

������

���Υ⥸�塼��ϡ�
netCDF�ե����ޥå�(http://www.unidata.ucar.edu/packages/netcdf/)
�ǥǡ�������Ϥ��뤿��Τ�ΤǤ���

----------------------------------------------------------------------

�����󥹥ȡ���

1. Makefile���Խ����ޥ���"CANS_LIB_DIR"���ͤ�Ŭ�������ꤹ�롣
2. make�ˤ�ꥳ��ѥ��롣

----------------------------------------------------------------------

������ˡ

CANS1D�ǥ����ȥ�ӥ塼������md_shktb�򻲹ͤˤ��롣

      call ncopn1d(idf,'out.cdf',idn,idi,ix) �ե�����򥪡��ץ󡢼��������
      call ncputattc(idf,'comment','cans1d')    �����Ȥ���ϡ�ʸ�����
      call ncdefss(idf,'t',idt,idn)       ��������ϥǡ���������ʥ����顼��ư��
      call ncdefs1(idf,'ro',idro,idn,idi)  ��������ϥǡ����������������ư��
      call ncdefs1(idf,'pr',idpr,idn,idi)  ��������ϥǡ����������������ư��
      call ncdefs1(idf,'vx',idvx,idn,idi)  ��������ϥǡ����������������ư��
      .......... 
      call ncputosi(idf,'ix',ix)          �����ϥǡ�������ϡʥ����顼������
      call ncputos(idf,'gm',gm)           �����ϥǡ�������ϡʥ����顼��ư��
      call ncputo1(idf,idi,'x',x,ix)      �����ϥǡ�������ϡ�������ư��
      .......... 
      call ncputss(idf,idt,nd,time)       ������ǡ�������ϡʥ����顼��ư��
      call ncputs1(idf,idro,nd,ro,ix)     ������ǡ�������ϡ�������ư��
      call ncputs1(idf,idpr,nd,pr,ix)     ������ǡ�������ϡ�������ư��
      call ncputs1(idf,idvx,nd,vx,ix)     ������ǡ�������ϡ�������ư��

      ..........
      ..........
      mstatus=nf_close(idf) �ǡ����ե�������Ĥ���

----------------------------------------------------------------------

���ƥ��֥롼����Ȱ�����������

ncopn1d(idf,file,idn,idi,ix)
ncopn2d(idf,file,idn,idi,ix,idj,jx)
ncopn3d(idf,file,idn,idi,ix,idj,jx,idk,kx)
  ��Ū
   �ǡ����ե�����򳫤��롣�ǡ����μ�����������롣
  ����
   idf    �ե����������ֹ�
   idn    �����ѿ����ۤˤ�������ʤ��ˤ������ֹ�
   idi,idj,idk    �����ѿ� ix,jx,kx�������ֹ�
  ����
   file  �ե�����̾
   ix,jx,kx    �ǡ������礭��

-------------------------------------------
ncdefss(idf,name,idda,idn)     ��ư�����������顼
ncdefs1(idf,name,idda,idn,idi) ��ư������1��������
ncdefs2(idf,name,idda,idn,idi,idj) ��ư������2��������
ncdefs3(idf,name,idda,idn,idi,idj,idk) ��ư������3��������
  ��Ū
   ��������ϥǡ��������
  ����
   idda   �ǡ����������ֹ�
  ����
   idf    �ե����������ֹ�
   idn    �����ѿ����ۤˤ�������ʤ��ˤ������ֹ�
   idi    �����ѿ� ix�������ֹ�
   idj    �����ѿ� jx�������ֹ�
   idk    �����ѿ� kx�������ֹ�
   name   �ѿ�̾��ʸ������

-------------------------------------------
ncputss(idf,idda,nd,da)       ��ư�����������顼
ncputs1(idf,idda,nd,da,ix)     ��ư����������1����
ncputs2(idf,idda,nd,da,ix,jx)     ��ư����������2����
ncputs3(idf,idda,nd,da,ix,jx,kx)     ��ư����������3����
 ��Ū
   ������ǡ��������
 ���Ϥʤ�
 ����
   idf �ե����������ֹ�
   idda   �ǡ����������ֹ�
   da �ǡ���
   ix,jx,kx    �ǡ������礭��
   nd �ǡ������Ǽ�����ֹ�

-------------------------------------------
ncputos(idf,name,da) ��ư�������������顼
ncputosi(idf,name,da) �����������顼
ncputo1(idf,idi,name,da,ix) ��ư��������1��������
ncputo2(idf,idi,idj,name,da,ix,jx) ��ư��������2��������
ncputo3(idf,idi,idj,idk,name,da,ix,jx,kx) ��ư��������3��������
  ��Ū
    �����ϥǡ��������
  ���Ϥʤ�
  ����
   idf    �ե����������ֹ�
   da  �ǡ���
   name   �ѿ�̾��ʸ������
   ix,jx,kx    �ǡ������礭��
   idi    �����ѿ� ix�������ֹ�
   idj    �����ѿ� jx�������ֹ�
   idk    �����ѿ� kx�������ֹ�
----------------------------------------------------------------------

