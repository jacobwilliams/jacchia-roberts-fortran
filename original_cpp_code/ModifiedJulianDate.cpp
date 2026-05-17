Real ModifiedJulianDate(YearNumber year, MonthOfYear month, DayOfMonth day,
                        Integer hour, Integer minute, Real second,
                        Real refEpochJD)
{
    // 5/23/06 - commented out
    // Need to subtract out Julian date offset before adding the fraction
    // of a day term to gain significant digits after the decimal place.
    // The computation is copied from the Julian Date method.

   Integer computeYearMon = ( 7*(year + (Integer)((month + 9)/12)) )/4;
   Integer computeMonth = (275 * month)/9;
   Real fractionalDay = ((second/60.0 + minute)/60.0 + hour)/24.0;

   Real ModJulianDay = 367*year - computeYearMon + computeMonth + day +
               1721013.5 -  refEpochJD;
   Real modJulianDate = ModJulianDay  + fractionalDay;

   return modJulianDate;
}
