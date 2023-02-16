
/*
 * C Subroutine that may be called from the basis package shl.
 *
 * 
 * */

int bar(argc,argv)
int argc;	/* 	number of arguments  */
void *argv[];	/* 	array of pointers to the arguments  */
{

	int i;
	float *a;	/* 	1st argument   */
	int *len;	/* 	2nd argument  */

	if(argc < 2){
		printf("usage: status = bar(<shared library>,\"bar\",<array>,<len>)\n");
		return 1;
	}
	a = (float *)argv[0];
	len = (int *)argv[1];

	for(i = 0; i < *len; ++i)a[i] = (float) i;
	return 0;

}
