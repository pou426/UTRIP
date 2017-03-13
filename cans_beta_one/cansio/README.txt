

CANS nc モジュール

       ver.0  2001.9.7

----------------------------------------------------------------------

◯概要

このモジュールは、
netCDFフォーマット(http://www.unidata.ucar.edu/packages/netcdf/)
でデータを出力するためのものです。

----------------------------------------------------------------------

◯インストール

1. Makefileを編集。マクロ"CANS_LIB_DIR"の値を適当に設定する。
2. makeによりコンパイル。

----------------------------------------------------------------------

◯使用法

CANS1Dディストリビューションのmd_shktbを参考にする。

      call ncopn1d(idf,'out.cdf',idn,idi,ix) ファイルをオープン、次元を定義
      call ncputattc(idf,'comment','cans1d')    コメントを出力（文字列）
      call ncdefss(idf,'t',idt,idn)       時系列出力データを定義（スカラー浮動）
      call ncdefs1(idf,'ro',idro,idn,idi)  時系列出力データを定義（配列浮動）
      call ncdefs1(idf,'pr',idpr,idn,idi)  時系列出力データを定義（配列浮動）
      call ncdefs1(idf,'vx',idvx,idn,idi)  時系列出力データを定義（配列浮動）
      .......... 
      call ncputosi(idf,'ix',ix)          一回出力データを出力（スカラー整数）
      call ncputos(idf,'gm',gm)           一回出力データを出力（スカラー浮動）
      call ncputo1(idf,idi,'x',x,ix)      一回出力データを出力（配列浮動）
      .......... 
      call ncputss(idf,idt,nd,time)       時系列データを出力（スカラー浮動）
      call ncputs1(idf,idro,nd,ro,ix)     時系列データを出力（配列浮動）
      call ncputs1(idf,idpr,nd,pr,ix)     時系列データを出力（配列浮動）
      call ncputs1(idf,idvx,nd,vx,ix)     時系列データを出力（配列浮動）

      ..........
      ..........
      mstatus=nf_close(idf) データファイルを閉じる

----------------------------------------------------------------------

◯各サブルーチンと引数の説明。

ncopn1d(idf,file,idn,idi,ix)
ncopn2d(idf,file,idn,idi,ix,idj,jx)
ncopn3d(idf,file,idn,idi,ix,idj,jx,idk,kx)
  目的
   データファイルを開ける。データの次元を定義する。
  出力
   idf    ファイル装置番号
   idn    時間変数（陽には定義しない）の装置番号
   idi,idj,idk    空間変数 ix,jx,kxの装置番号
  入力
   file  ファイル名
   ix,jx,kx    データの大きさ

-------------------------------------------
ncdefss(idf,name,idda,idn)     浮動小数点スカラー
ncdefs1(idf,name,idda,idn,idi) 浮動小数点1次元配列
ncdefs2(idf,name,idda,idn,idi,idj) 浮動小数点2次元配列
ncdefs3(idf,name,idda,idn,idi,idj,idk) 浮動小数点3次元配列
  目的
   時系列出力データを定義
  出力
   idda   データの装置番号
  入力
   idf    ファイル装置番号
   idn    時間変数（陽には定義しない）の装置番号
   idi    空間変数 ixの装置番号
   idj    空間変数 jxの装置番号
   idk    空間変数 kxの装置番号
   name   変数名（文字型）

-------------------------------------------
ncputss(idf,idda,nd,da)       浮動小数点スカラー
ncputs1(idf,idda,nd,da,ix)     浮動小数点配列1次元
ncputs2(idf,idda,nd,da,ix,jx)     浮動小数点配列2次元
ncputs3(idf,idda,nd,da,ix,jx,kx)     浮動小数点配列3次元
 目的
   時系列データを出力
 出力なし
 入力
   idf ファイル装置番号
   idda   データの装置番号
   da データ
   ix,jx,kx    データの大きさ
   nd データを格納する番号

-------------------------------------------
ncputos(idf,name,da) 浮動小数点型スカラー
ncputosi(idf,name,da) 整数型スカラー
ncputo1(idf,idi,name,da,ix) 浮動小数点型1次元配列
ncputo2(idf,idi,idj,name,da,ix,jx) 浮動小数点型2次元配列
ncputo3(idf,idi,idj,idk,name,da,ix,jx,kx) 浮動小数点型3次元配列
  目的
    一回出力データを出力
  出力なし
  入力
   idf    ファイル装置番号
   da  データ
   name   変数名（文字型）
   ix,jx,kx    データの大きさ
   idi    空間変数 ixの装置番号
   idj    空間変数 jxの装置番号
   idk    空間変数 kxの装置番号
----------------------------------------------------------------------

