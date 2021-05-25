# INFO MissingValuesMeteo
Code to imput missing temperatures data, 3 ways to generate the missing data:

* Option 1: na.interp (a linear way)
* Option 2: using the package "imputeTS" based on the Kalman filter
* Option 3: ARIMA

here we can compare the results of the 3 options above:

| Option      |  2005-04-08 21:00:00  |  2005-04-11 00:00:00  |  2005-04-11 03:00:00 |
| ----------- | ----------- | ----------- | ----------- |
| na.interp   | 10.5       | 11       | 10       |
| imputeTS    | 10.203376  | 9.926740 | 8.567582 |
| ARIMA       | 12.492506  | 8.909408 | 8.097138 |
