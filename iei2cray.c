/* convert 32-bit IEEE integer to 64-bit Cray integer */
void IEI2CRAY(myieee, mycray, nword)
int *myieee, *mycray, *nword;
{
	int iword = *nword, sign;
	int *ieeeword, *crayword;
	while (iword > 0) {
		iword--;
		/* pointer to IEEE word */
		ieeeword = myieee + iword / 2;
		/* pointer to Cray word */
		crayword = mycray + iword;
		/* copy IEEE to Cray */
		*crayword = *ieeeword;
		/* if second IEEE word shift 4 bytes left */
		if (iword % 2 == 1) *crayword = *crayword << 32;
		/* shift 4 bytes right */
		*crayword = *crayword >> 32;
		/* find out the sign from bit 31 */
		sign = *crayword >> 31;
		/* if negative replace leftmost 32 0's with 1's */
		if (sign == 1) *crayword = *crayword | 0xFFFFFFFF00000000;
	}
	return;
}
