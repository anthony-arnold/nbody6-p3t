#define NFIELD 3
#define NLINMAX 30

int nlin[NFIELD]={NLINMAX, 26, 9};
char *namef[NFIELD] = {"Sarajedini","Simioni", "Sohn"};

double ragrid[NFIELD][NLINMAX]={{47.99631920, 47.99824769, 47.99737310, 48.04364810, 48.04477211, 48.04541299, 48.04450820, 48.08940670, 48.09302296, 48.13461620, 48.13265643, 48.13355820, 48.13159936, 48.13250120, 48.13054156, 48.13144420, 48.12948534, 48.13038720, 48.08554010, 48.08497350, 48.08316032, 48.08403200, 48.03781110, 47.99316020, 47.99508787, 47.99421320, 47.99614082, 47.99526620, 47.99719312, 47.99631920}, {48.31460200, 48.24575600, 48.16926300, 48.23809400,  48.23723400,  48.23393500, 48.13051900, 48.12699100, 48.05059000, 47.98138100, 47.98323400, 47.88062400, 47.88630500, 47.82164000, 47.89602500, 47.96705200, 47.89870000, 47.98876300, 47.98323400, 47.98138100, 48.05776500, 48.12699100, 48.13051900, 48.13397000, 48.23723400, 48.31460200},
 {48.1296158,  48.0517197,  48.0699158,  48.0699196,  48.0699501,  48.1486816,  48.1486435,  48.1486435,  48.1296158}};

double decgrid[NFIELD][NLINMAX]={{-55.23334790, -55.23380940, -55.23480030, -55.24586190, -55.24458663, -55.24473639, -55.24572960, -55.25621060, -55.25223540, -55.20644880, -55.20599205, -55.20499810, -55.20454153, -55.20354750, -55.20309071, -55.20209580, -55.20163915, -55.20064510, -55.19017910, -55.19080281, -55.19036995, -55.18937820, -55.17833180, -55.22899190, -55.22945331, -55.23044430, -55.23090566, -55.23189660, -55.23235776, -55.23334790}, {-55.22560900, -55.18226400, -55.22174900, -55.26513600, -55.26245500, -55.32091100, -55.31897500, -55.31796200, -55.35713000, -55.31336900, -55.31089100, -55.30771600, -55.24894400, -55.21041800, -55.16967700, -55.21184700, -55.24935000, -55.25211400, -55.31089100, -55.31336900, -55.27424400, -55.31796200, -55.31897500, -55.26052200, -55.26245500,  -55.22560900},
 { -55.2875252, -55.2948341, -55.3394356, -55.3394356, -55.3395157, -55.3315964, -55.3315125, -55.3315125, -55.2875252}};

#define FITFILE "../ngc1261_bestfit"
