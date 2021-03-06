インストールと実行

CANS2D ディストリビューション

       ver.0  2001.9.7 

----------------------------------------------------------------------

◯概要

このディストリビューションは、流体・MHDなどの数値計算のための
Fortranモジュール群です。計算例もモデルパッケージとしてはいっています。

このディストリビューションには2次元計算用のコードが含まれています。

----------------------------------------------------------------------
◯動作環境

FreeBSD4.2
Linux 2.2
IRIX64 6.5
SunOS5.7 'make "FC=f90"'
HP-UX 11.0 'make "FC=f90"'

SUPER-UX 11.1 (SX5) 'make "FC=f90"'
UXP/V (VPP) 'make "FC=frt"'


----------------------------------------------------------------------
◯インストールと実行

0. netCDFをインストールする。

1. makeによりモジュール群をコンパイル。
   ひとつ上のディレクトリにlibcans2d.aというファイルができる。

2. モデルパッケージで計算を試す。ディレクトリmd_*のどれかの下で
   makeをおこなう。計算が実行されて結果がファイルout.cdfに出力される。

3. IDLでデータを読込んで（rddt.proを使う）可視化する。CANS
   ディストリビューションのディレクトリidl/をIDL_PATHに加えるのを
   忘れずに。

----------------------------------------------------------------------

◯各ディレクトリの説明

CANS2D モジュール

  common/ 計算で使うさまざまな共通ルーチンを集めたモジュール
  bc/     境界条件を定義するためのモジュール

  hdmlw/  流体力学方程式・MHD方程式を
             改良Lax-Wendroff＋人工粘性法で解くためのモジュール
  hdroe/  流体力学方程式・MHD方程式を
             Roe＋TVD法で解くためのモジュール

  cndsor/ 熱伝導を陰解法（時間精度1次：行列反転はRed Black SOR法）で解く
              ためのモジュール
  cndbicg/ 熱伝導を陰解法（時間精度1次：行列反転はbiCG-stab法）で解く
              ためのモジュール
  htcl/   放射冷却・静的加熱を陽解法で解くためのモジュール
  selfg/  自己重力を解くためのモジュール

CANS2D モデルパッケージ

  md_shktb/   流体衝撃波菅問題   [流体]
  md_sedov/   球対称点源爆発（Sedov解）問題  [非一様断面積、流体]
  md_cndtb/   単純熱伝導問題  [熱伝導]
  md_cndsp/   球対称単純熱伝導問題  [非一様断面積、熱伝導]
  md_thinst/  熱不安定問題   [冷却加熱]
  md_flare/   フレアループ問題  [非一様断面積、流体、熱伝導、冷却加熱]

  md_mhdshktb/ MHD衝撃波菅問題  [MHD]
  md_spicule/ 磁力管Alfven波伝播（スピキュール生成）問題 [非一様断面積、MHD]

----------------------------------------------------------------------
