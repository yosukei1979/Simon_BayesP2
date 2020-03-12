/*******************************************************************************
    FILE NAME    : bayesP2sim.sas
    TITLE   : Bayesian phase II trial design by posterior probability
    PRODUCT : SAS R9.4
    AUTHOR  : Y.Inaba
    DATE    : 2020/3/12
*******************************************************************************/
proc iml;
/*******************************************************************************
    例数設計の確認
*******************************************************************************/
**サンプルサイズ;
n_max=10;
*y=8;
***標準治療の有効確率;
p_s=0.4;

***p_sの従う事前分布のパラメータ(Beta(alpha_s, beta_s));
alpha_s=14;
beta_s =21;
*CALL RANDGEN (p_e, "Beta" , alpha_s, beta_s);
p_e=0.8;
***試験治療の有効確率;


***p_eの従う事前分布のパラメータ(Beta(alpha_e, beta_e));
alpha_e=1;
beta_e =1;

***反応率の最低限許容可能な増分;
delta=0;
delta2=1-delta;

/*result=j(1,1,0);*/
 
**非積分関数の定義;
start func(p) global(alpha_s, beta_s, alpha_e, beta_e, delta, n, y);
	F=CDF("BETA", p+delta, alpha_e+y, beta_e+n-y);
	f2=PDF("BETA", p, alpha_s, beta_s);
	res_=(1-F)*f2;
return(res_);
finish;
**積分範囲;
a={0} || delta2;

**標本サイズ検討結果の確認;
/*do n=6 to n_max;*/
/*	do y=0 to n;*/
/*		call quad(res,"func",a);*/
/*		_result=n||y||res;*/
/*		result=result//_result;*/
/*		free _result;*/
/*	end;*/
/*end;*/
/**/
/*print result;*/
/*free result;*/


/*******************************************************************************
    帰無仮説が正しい時に試験成功する確率を算出（Type I error）
*******************************************************************************/
***--------------------
   n=10固定の場合
***--------------------;

**Trial数;
t=10000;
n=10;
**帰無仮説が真;
p_e=0.4;

/*x=j(t,3,0);*/
/*do i=1 to t;*/
/*	y = rand("Binomial", p_e, n);*/
/*	call quad(res2,"func",a);*/
/*	x[i, 1] = y;*/
/*	x[i, 2] = res2;*/
/*	x[i, 3] = res2>0.95;*/
/*end;*/
/**/
/*res3=loc(x[,3]);*/
/**/
/*res4=ncol(res3)/t;*/
/*print t, res4;*/
/*free t,n;*/
***--------------------
   n=7,8,9,10どこかで確率を超えたら終了するデザイン
***--------------------;
t=10000;


x=j(t,9,0);

do i=1 to t;

	**n=7;
	y = rand("Binomial", p_e, 7);
	call quad(res2,"func",a);
	x[i, 1] = y;
	x[i, 2] = res2;

	**n=8;
	y = y + rand("Binomial", p_e, 1);
	call quad(res2,"func",a);
	x[i, 3] = y;
	x[i, 4] = res2;

	**n=9;
	y = y + rand("Binomial", p_e, 1);
	call quad(res2,"func",a);
	x[i, 5] = y;
	x[i, 6] = res2;

	**n=10;
	y = y + rand("Binomial", p_e, 1);
	call quad(res2,"func",a);
	x[i, 7] = y;
	x[i, 8] = res2;

	x[i, 9] = x[i, 2]>0.99 | x[i, 4]>0.99 | x[i, 6]>0.99 | x[i, 8]>0.95 ;
	
end;

res3=loc(x[,9]);

res4=ncol(res3)/t;
print t, res4;
*free t n;

quit;
