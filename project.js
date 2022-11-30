var canvas;
var gl;

var eye = [1, 1, 1];
var at = vec3(0.0, 0.0, 0.0);
var up = vec3(0.0, 1.0, 0.0);
var rotate_star = false;

var numtop  = 48;

var pointsArray = [];
var normalsArray = [];
var texCoordsArray = [];
var textures = [];

var tableIndex;
var tableNum;

var starIndex;
var starNum;

var firstIndex
var firstNum
var secondIndex;
var secondNum;
var topIndex;
var topNum;

var gingerbreadIndex;
var gingerbreadNum;

var treeIndex;
var treeNum;

var lightPosition = vec4(3.0, -3.0, -2.0, 1.0 );
var lightAmbient = vec4(1.0, 1.0, 1.0, 1.0 );
var lightDiffuse = vec4( 1.0, 1.0, 1.0, 1.0 );
var lightSpecular = vec4( 1.0, 1.0, 1.0, 1.0 );

var ctm;
var ambientColor, diffuseColor, specularColor;
var modelView, projectionMatrix;
var modelViewLoc, projectionMatrixLoc; //currently just used in Bev
var program;
//404
var thetaLoc;

var AllInfo = {

    // Camera pan control variables.
    zoomFactor : 2,
    translateX : 0,
    translateY : 2,

    // Camera rotate control variables.
    phi : .8,
    theta : 0.5,
    radius : 1,
    dr : 2.0 * Math.PI/180.0,

    // Mouse control variables
    mouseDownRight : false,
    mouseDownLeft : false,

    mousePosOnClickX : 0,
    mousePosOnClickY : 0
};


var base = [ //Used in making the table
	vec4(-4, -1, 1, 1), //A 0
	vec4(-4, -1, 0, 1), //B 1
	vec4(-2, -1, 0, 1), //C 2
	vec4(-2, -1, 1, 1), //D 3
]

var treePoints = [
  //top
	[0,    .05, 0.0],
	[.034, .115, 0.0],
	[.045, .135, 0.0],
	[.065, .170, 0.0],
  //end of top
  //Second set start
  [.035, .176, 0.0],
  [.050, .190, 0.0],
	[.080, .225, 0.0],
  //end of second
  //beginning of third
	[.055, .230, 0.0],
	[.065, .240, 0.0],
	[.080, .257, 0.0],
	[.090, .270, 0.0],
	[.100, .280, 0.0],
  //end of third
  //start of fourth
	[.075, .285, 0.0],
	[.090, .305, 0.0],
	[.104, .320, 0.0],
	[.130, .350, 0.0],
  //end of third
  //start of fourth
	[.105, .355, 0.0],
	[.120, .380, 0.0],
	[.160, .433, 0.0],
  //end of fourth
  //start of fifth
	[.135, .435, 0.0],
	[.145, .455, 0.0],
	[.155, .468, 0.0],
	[.185, .512, 0.0],
	[.030, .526, 0.0],
	[0, .525, 0.0],
];

function scale4(a, b, c) {
        var result = mat4();
        result[0][0] = a;
        result[1][1] = b;
        result[2][2] = c;
        return result;
}

function MakeChristmasTree() {
	treeIndex = pointsArray.length;
	//Setup initial points matrix
	let top = [];
	for (let i = 0; i<25; i++)
	{
		top.push(vec4(treePoints[i][0], treePoints[i][1],
                                   treePoints[i][2], 1));
	}

	var r;
    var t=Math.PI/12;

        // sweep the original curve another "angle" degree
	for (let j = 0; j < 24; j++)
	{
                var angle = (j+1)*t;

                // for each sweeping step, generate 25 new points corresponding to the original points
		for(let i = 0; i < 25 ; i++ )
		{
		        r = top[i][0];
                        top.push(vec4(r*Math.cos(angle), top[i][1], -r*Math.sin(angle), 1));
		}
	}

       var N=25;

       for (let i=0; i<24; i++) // slices
       {
           for (let j=0; j<24; j++)  // layers
           {
						 //HELP//
				quad(i*N+j, (i+1)*N+j, (i+1)*N+(j+1), i*N+(j+1), top);
           }
       }
		treeNum = pointsArray.length - treeIndex;
}


function ExtrudedShape(top, height) {
  let bot = [];
  let N = top.length;
  let both = [];

	for (let i = 0; i < top.length; i++) {
		bot[i] = top[i].slice();
		both.push(top[i].slice());
	}

	//Offset the bottom points by height
	for (let i = 0; i < bot.length; i++) {
		bot[i][1] += height; //height;
		both.push(bot[i].slice());
	}
	for (let i = 0; i < N; i++) {
		quad(i, i + N, (i + 1 + N) % both.length, (i + 1) % N, both);
	}

	//Create the top and bottom polygons
	let indices = [];
	for(let i = 0; i <= top.length; i++) {
		indices.push(i % top.length);
	}
	polygon(indices, top);
	polygon(indices, bot);
}


function MakeGingerbreadMan() {
	gingerbreadIndex = pointsArray.length;
    let height=.1;
    let radius=.55;
    let num=30;
    var alpha= 2 * Math.PI/num ;

    let top = [vec4(0, 0, 0, 1)];

    for (let i=num; i>=0; i--)
    {
        top.push(vec4(radius*Math.cos(i*alpha), 0, radius*Math.sin(i*alpha), 1));
    }
    top.push(vec4(1.3, 0, 0, 1));
    top.push(vec4(1.3, 0, .15, 1));
    top.push(vec4(.5, 0, .2, 1));
    top.push(vec4(1, 0, .3, 1));
    top.push(vec4(1.8, 0, .5, 1));
    top.push(vec4(1.8, 0, .3, 1));
    top.push(vec4(1.3, 0, .4, 1));
    top.push(vec4(1.3, 0, .5, 1));
    top.push(vec4(1.8, 0, .6, 1));
    top.push(vec4(1.8, 0, .7, 1));
    top.push(vec4(1.8, 0, .8, 1));
    top.push(vec4(.6, 0, .4, 1));
    top.push(vec4(.6, 0, .45, 1));
    top.push(vec4(1.2, 0, .88, 1));
    top.push(vec4(1.2, 0, .65, 1));

    ExtrudedShape(top, height);
  	gingerbreadNum = pointsArray.length - gingerbreadIndex;
}


var floorIndex; var floorNum;
var wallIndex; var wallNum;
function MakeRoom() {
	floorIndex = pointsArray.length;
	let floor_corners = [ //Used in making the table
		vec4(-4, -1, 4, 1), //A 0
		vec4(-4, -1, -4, 1), //B 1
		vec4(4, -1, -4, 1), //C 2
		vec4(4, -1, 4, 1), //D 3
	]
	MakePrism(floor_corners[0], floor_corners[1], floor_corners[2], floor_corners[3], .1);

	floorNum = pointsArray.length - floorIndex;
	wallIndex = pointsArray.length;
	let left = [
		vec4(-4, 5, 4, 1), //4 0
		vec4(-4, 5, -4, 1), //B 1
		vec4(-4.1, 5, -4, 1), //4 0
		vec4(-4.1, 5, 4, 1), //B 1
	];
	let back = [
		vec4(4, 5, -4, 1), //B 1
		vec4(-4, 5, -4, 1), //C 2
		vec4(-4, 5, -4.1, 1), //B 1
		vec4(4, 5, -4.1, 1), //C 2
	];
	let right = [
		vec4(4, 5, 4, 1),
		vec4(4, 5, -4, 1),
		vec4(4.1, 5, -4, 1),
		vec4(4.1, 5, 4, 1),
	];

	MakePrism(left[0], left[1], left[2], left[3], -6);
	MakePrism(back[0], back[1], back[2], back[3], -6);
	//MakePrism(right[0], right[1], right[2], right[3], -6);

	wallNum = pointsArray.length - wallIndex;
}


function MakeOrnament() {
	ornamentIndex = pointsArray.length;
	var va = vec4(0.0, 0.0, -1.0,1);
	var vb = vec4(0.0, 0.942809, 0.33333, 1);
	var vc = vec4(-0.816497, -0.471405, 0.333333, 1);
	var vd = vec4(0.816497, -0.471405, 0.333333,1);
	var numTimesToSubdivide = 6;
	tetrahedron(va, vb, vc, vd, numTimesToSubdivide);

	ornamentNum = pointsArray.length - ornamentIndex;

	topIndex = pointsArray.length;

	topNum = pointsArray.length - topIndex;

}

//----------------------------- BEV STUFF------------------------------------------//

function HalfCircle() { //USED IN ORNAMENT
    var height=.25;
    var radius=.8;
    var num=30;
    var alpha= 2 * Math.PI/num;

    var top = [vec4(0, 0, 0, 1)];
    for (var i=num; i>=0; i--)
    {
        top.push(vec4(radius*Math.cos(i*alpha), 0, radius*Math.sin(i*alpha), 1));
    }

    N=N_Circle=top.length;

    // add the second set of points
    for (var i=0; i<N; i++)
    {
        top.push(vec4(top[i][0], top[i][1]+height, top[i][2], 1));
    }
    ExtrudedShape(top);
}


function polygon(indices, shape) { //makes a polygon
    // for indices=[a, b, c, d, e, f, ...]
    var M=indices.length;
    var normal=Newell(indices, shape);

    var prev=1;
    var next=2;
    for (var i=0; i<M-2; i++) {
        pointsArray.push(shape[indices[0]]);
        normalsArray.push(normal);

        pointsArray.push(shape[indices[prev]]);
        normalsArray.push(normal);

        pointsArray.push(shape[indices[next]]);
        normalsArray.push(normal);

		texCoordsArray.push(vec2(0, 0));
		texCoordsArray.push(vec2(0, 1));
		texCoordsArray.push(vec2(1, 1));

        prev=next;
        next=next+1;
    }
}

function Newell(indices, shape) { //idk what this is lol
   var L=indices.length;
   var x=0, y=0, z=0;
   var index, nextIndex;

   for (var i=0; i<L; i++)
   {
       index=indices[i];
       nextIndex = indices[(i+1)%L];

       x += (shape[index][1] - shape[nextIndex][1])*
            (shape[index][2] + shape[nextIndex][2]);
       y += (shape[index][2] - shape[nextIndex][2])*
            (shape[index][0] + shape[nextIndex][0]);
       z += (shape[index][0] - shape[nextIndex][0])*
            (shape[index][1] + shape[nextIndex][1]);
   }

   return (normalize(vec3(x, y, z)));
}

function triangle(a, b, c) {
     var t1 = subtract(b, a);
     var t2 = subtract(c, a);
     var normal = normalize(cross(t1, t2));

     normalsArray.push(normal);
     normalsArray.push(normal);
     normalsArray.push(normal);

     pointsArray.push(a);
     pointsArray.push(b);
     pointsArray.push(c);

	texCoordsArray.push(vec2(0, 0));
    texCoordsArray.push(vec2(0, 1));
    texCoordsArray.push(vec2(1, 1));
}

function divideTriangle(a, b, c, count) { //this code turns a triangle into a sphere
    if ( count > 0 ) {
        var ab = mix( a, b, 0.5);
        var ac = mix( a, c, 0.5);
        var bc = mix( b, c, 0.5);

        ab = normalize(ab, true);
        ac = normalize(ac, true);
        bc = normalize(bc, true);

        divideTriangle( a, ab, ac, count - 1 );
        divideTriangle( ab, b, bc, count - 1 );
        divideTriangle( bc, c, ac, count - 1 );
        divideTriangle( ab, bc, ac, count - 1 );
    }
    else {
        triangle( a, b, c );
    }
}

function tetrahedron(a, b, c, d, n) { //make a tetrahedron. n is the num of times to subdivide
    divideTriangle(a, b, c, n);
    divideTriangle(d, c, b, n);
    divideTriangle(a, d, b, n);
    divideTriangle(a, c, d, n);
}

//----------------------------- END -----------------------------------------------------------------//

//---------------------------------- ERIC OBJECTS ---------------------------------------------------//
function MakeTable()
{
	tableIndex = pointsArray.length;
	MakePrism(base[0], base[1], base[2], base[3], .1);
	//Now the legs
	var A, B, C, D;
	var w = .1;
	var h = .9;
	//Front Left
	A = base[0];
	B = vec4(A[0], A[1], A[2] - w, 1.0);
	C = vec4(B[0] + w, B[1], B[2]);
	D = vec4(C[0], C[1], C[2] + w);
	//MakePrism(A, B, C, D, h);
	MakeCylinder(Average(A, B, C, D), w/2, h);
	//Back Left
	B = base[1];
	A = vec4(B[0], B[1], B[2] + w, 1.0);
	C = vec4(B[0] + w, B[1], B[2]);
	D = vec4(C[0], C[1], C[2] + w);
	//MakePrism(A, B, C, D, h);
	MakeCylinder(Average(A, B, C, D), w/2, h);
	//Back Right
	C = base[2];
	D = vec4(C[0], C[1], C[2] + w);
	A = vec4(D[0] - w, D[1], D[2]);
	B = vec4(A[0], A[1], A[2] - w, 1.0);
	//MakePrism(A, B, C, D, h);
	MakeCylinder(Average(A, B, C, D), w/2, h);
	//Front Right
	D = base[3];
	A = vec4(D[0] - w, D[1], D[2]);
	B = vec4(A[0], A[1], A[2] - w, 1.0);
	C = vec4(B[0] + w, B[1], B[2]);
	//MakePrism(A, B, C, D, h);
	MakeCylinder(Average(A, B, C, D), w/2, h);
	tableNum = pointsArray.length - tableIndex;
}

var cupIndex;
var cupNum;
var milkIndex;
var milkNum;

function MakeMilk() {
	cupIndex = pointsArray.length;
	let top = [];
	let bot = [];
	let both_outer = [];
	let both_inner = [];
	let angle = 0;
	let radius = .15;
	let center = [-3, -.5, .5, 1];
	let height = .3;
	let num = 64;
	for(let j = 0; j <= num+1; j++ ){
		angle = j * (2.0*Math.PI/num);
		top.push(vec4(Math.cos(angle) * radius + center[0], center[1], Math.sin(angle) * radius + center[2])); //AT
		top.push(vec4(Math.cos(angle) * radius*.95 + center[0], center[1], Math.sin(angle) * radius*.95 + center[2])); //BT
		bot.push(vec4(Math.cos(angle) * radius + center[0], center[1]-height, Math.sin(angle) * radius + center[2])); //AB
		bot.push(vec4(Math.cos(angle) * radius*.95 + center[0], center[1]-height, Math.sin(angle) * radius*.95 + center[2])); //BB
		both_outer.push(top[top.length-1]);
		both_outer.push(bot[bot.length-1]);
		both_inner.push(top[top.length-2]);
		both_inner.push(bot[bot.length-2]);
		if (j > 1) {
			quad(top.length-4, top.length-3, top.length-1, top.length-2, top);
			quad(bot.length-4, bot.length-3, bot.length-1, bot.length-2, bot);
			//sides
			quad(both_outer.length-4, both_outer.length-3, both_outer.length-1, both_outer.length-2, both_outer);
			quad(both_inner.length-4, both_inner.length-3, both_inner.length-1, both_inner.length-2, both_inner);
		}
	}
	quad(3, 2, 0, 1, top);
	quad(3, 2, 0, 1, bot);
	center[1] -= height;
	MakeCylinder(center, radius * .95, height * .1, num);
	cupNum = pointsArray.length - cupIndex;

	center[1] += height * .7;
	milkIndex = pointsArray.length;
	MakeCylinder(center, radius * .949, height * .6, num);
	milkNum = pointsArray.length - milkIndex;
}

var plateIndex; var plateNum;
var cookieIndex; var cookieNum;
var chipIndex; var chipNum;
function MakeCookies() {
	//what have i done
	plateIndex = pointsArray.length; //make the plate
	let top_center = [0, 0, 0, 0];
	let top_radius = 3;
	let bot_radius = top_radius/2;
	let bot_center = [0, -top_radius/6, 0, 0];
	MakePlate(top_center, bot_center, top_radius, bot_radius);
	plateNum = pointsArray.length - plateIndex;


	cookieIndex = pointsArray.length; //make the cookie

	let init = 3;
	num = 6;
	top_radius = init;
	top_center = [0, 0-init, 0, 0];

	for (let i = 0; i < num; i++) {
		bot_center = top_center.slice();
		bot_radius = top_radius;
		top_radius = init * (1-(i*1/num));
		top_center = [(-(i%2)/(init*5/3)), bot_center[1], ((i%2)/(init*5/3)), 0];
		top_center[1] += init/num * (1-(i*1/num));
		MakeConeThing(top_center, bot_center, top_radius, bot_radius, Math.floor(32 - (16 * (i/num))));
	}

	cookieNum = pointsArray.length - cookieIndex;

	chipIndex = pointsArray.length; //make the chocolate chips
	let y = -1 * init;
	let chip_point = [];
	for (let i = 0; i < num; i++) {
		y += init/num * (1-(i*1/num));
		if (i > 0) {
			for (let j = i/3; j < num/2; j++) {
				chip_point[0] = (-(i%2)/(init*5/3)) + Math.cos(j * (2.0*Math.PI/(num/2))) * (init*.95 * (1-(Math.max(1, i)*1/num)));
				chip_point[1] = y;
				chip_point[2] = ((i%2)/(init*5/3)) + Math.sin(j * (2.0*Math.PI/(num/2))) * (init*.95 * (1-(Math.max(1, i)*1/num)));
				MakeCone(chip_point.slice(), init/15, init/10);
			}
		}
	}
	chipNum = pointsArray.length -  chipIndex;
}

var starPoints = [ //10 points
	vec4(0, -.2, 0, 1), //A
	vec4(.4*Math.cos(54*Math.PI/180), -.4*Math.sin(54*Math.PI/180), 0, 1), //B
	vec4(.2*Math.cos(18*Math.PI/180), -.2*Math.sin(18*Math.PI/180), 0, 1), //C
	vec4(.4*Math.cos(18*Math.PI/180), .4*Math.sin(18*Math.PI/180), 0, 1), //D
	vec4(.2*Math.cos(54*Math.PI/180), .2*Math.sin(54*Math.PI/180), 0, 1), //E
	vec4(0, .4, 0, 1), //F
	vec4(-.2*Math.cos(54*Math.PI/180), .2*Math.sin(54*Math.PI/180), 0, 1), //G
	vec4(-.4*Math.cos(18*Math.PI/180), .4*Math.sin(18*Math.PI/180), 0, 1), //H
	vec4(-.2*Math.cos(18*Math.PI/180), -.2*Math.sin(18*Math.PI/180), 0, 1), //I
	vec4(-.4*Math.cos(54*Math.PI/180), -.4*Math.sin(54*Math.PI/180), 0, 1), //J
	vec4(0, 0, .1, 1) //CENTER
];

function MakeStar() {
	starIndex = pointsArray.length;
	//Point, Next Point, Center, Point
	for (let i = 0; i < starPoints.length; i++) {
		starPoints[i][0] -= 3;
	}
	for (let j = 0; j < 2; j++) {
		if (j == 1) {
			starPoints[starPoints.length-1] = vec4(-3, 0, -.1, 1);
		}
		for (let i = 0; i < starPoints.length; i++) {
			if (j == 0) {
				quad(i, starPoints.length-1, (i+1) % (starPoints.length-1), i, starPoints);
			}
			else {
				quad(i, (i+1) % (starPoints.length-1), starPoints.length-1, i, starPoints);
			}
		}
	}
	starNum = pointsArray.length - starIndex;
}

var present = [
	vec4( -0.5, -0.5,  0.5, 1.0 ), //A (0)
	vec4( -0.5,  0.5,  0.5, 1.0 ), //B (1)
	vec4( 0.5,  0.5,  0.5, 1.0 ), //C (2)
	vec4( 0.5, -0.5,  0.5, 1.0 ), //D (3)
	vec4( -0.5, -0.5, -0.5, 1.0 ), //E (4)
	vec4( -0.5,  0.5, -0.5, 1.0 ), //F (5)
	vec4( 0.5,  0.5, -0.5, 1.0 ), //G (6)
	vec4( 0.5, -0.5, -0.5, 1.0 ) //H (7)
];

function MakePresent() {
	presentIndex = pointsArray.length;
	quad( 1, 0, 3, 2, present);
	quad( 2, 3, 7, 6, present);
	quad( 3, 0, 4, 7, present);
	quad( 6, 5, 1, 2, present);
	quad( 4, 5, 6, 7, present);
	quad( 5, 4, 0, 1, present);
	presentNum = pointsArray.length - presentIndex;
}

function Average(A, B, C, D, xOff = 0, yOff = 0, zOff = 0) {
	var point = vec4(0, 0, 0, 1);
	point[0] = (A[0] + B[0] + C[0] + D[0])/4 + xOff;
	point[1] = (A[1] + B[1] + C[1] + D[1])/4 + yOff;
	point[2] = (A[2] + B[2] + C[2] + D[2])/4 + zOff;
	return point;
}

function MakePlate(top_center, bot_center, top_radius, bot_radius, num = 32){
	top_center[3] = 1;
	bot_center[3] = 1;
	let mid_center = bot_center.slice();
	mid_center[1] += (top_center[1] - bot_center[1])/2;
	var top = [top_center];
	var bot = [bot_center];
	var mid = [mid_center];
	var corners = []; //connects plate ring to bottom circle, ie bottom of plate
	var sec_corners = []; //connects plate ring to to p circle, ie top side of plate
	let angle = 0;

	for(let j = 0; j <= num; j++ ){
		angle = j * (2.0*Math.PI/num);
		top.push(vec4(Math.cos(angle) * top_radius + top_center[0], top_center[1], Math.sin(angle) * top_radius + top_center[2]));
		bot.push(vec4(Math.cos(angle) * bot_radius + bot_center[0], bot_center[1], Math.sin(angle) * bot_radius + bot_center[2]));
		mid.push(vec4(Math.cos(angle) * bot_radius + mid_center[0], mid_center[1], Math.sin(angle) * bot_radius + mid_center[2]));
		if (j > 0) {
			triangle(bot[0], bot[bot.length-1], bot[bot.length-2]);
			triangle(mid[mid.length-2], mid[mid.length-1], mid[0]);
			corners.push(top[top.length-1]);
			corners.push(bot[bot.length-1]);
			sec_corners.push(top[top.length-1]);
			sec_corners.push(mid[mid.length-1]);
		}
	}
	for (let i = 0; i < num*2; i++) {
		quad((0+i)%(2*num), (1+i)%(2*num), (3+i)%(2*num), (2+i)%(2*num), corners);
	}
	for (let i = 0; i < num*2; i++) {
		quad((0+i)%(2*num), (1+i)%(2*num), (3+i)%(2*num), (2+i)%(2*num), sec_corners);
	}
}

function MakeCone(center, radius, height, num = 32) {
	center[3] = 1;
	let bot = [center];
	let top = center.slice();
	top[1] += height;
	let angle = 0;

	for(let j = 0; j <= num; j++ ){
		angle = j * (2.0*Math.PI/num);
		bot.push(vec4(Math.cos(angle) * radius + center[0], center[1], Math.sin(angle) * radius + center[2]));
		if (j > 0) {
			triangle(bot[bot.length-1], bot[bot.length-2], bot[0]);
		}
	}
	for (let i = 0; i < bot.length; i++) {
		triangle(top, bot[i%bot.length], bot[(i+1)%bot.length]);
	}
}

function MakeConeThing(top_center, bot_center, top_radius, bot_radius, num = 32){
	top_center[3] = 1;
	bot_center[3] = 1;
	var top = [top_center];
	var bot = [bot_center];
	var corners = [];
	let angle = 0;

	for(let j = 0; j <= num; j++ ){
		angle = j * (2.0*Math.PI/num);
		top.push(vec4(Math.cos(angle) * top_radius + top_center[0], top_center[1], Math.sin(angle) * top_radius + top_center[2]));
		bot.push(vec4(Math.cos(angle) * bot_radius + bot_center[0], bot_center[1], Math.sin(angle) * bot_radius + bot_center[2]));
		if (j > 0) {
			triangle(top[top.length-1], top[0], top[top.length-2]);
			triangle(bot[0], bot[bot.length-1], bot[bot.length-2]);
			corners.push(top[top.length-1]);
			corners.push(bot[bot.length-1]);
		}
	}
	for (let i = 0; i < num*2; i++) {
		quad((0+i)%(2*num), (1+i)%(2*num), (3+i)%(2*num), (2+i)%(2*num), corners);
	}
}

function MakeCylinder(center, radius, height, num = 32) {
	var top_center = [center[0], center[1], center[2], 1]
	var top = [top_center];
	var bot_center = [center[0], center[1]-height, center[2], 1]
	var bot = [bot_center];
	var corners = [];
	let angle = 0;

	for(let j = 0; j <= num; j++ ){
		angle = j * (2.0*Math.PI/num);
		top.push(vec4(Math.cos(angle) * radius + center[0], center[1], Math.sin(angle) * radius + center[2]));
		bot.push(vec4(Math.cos(angle) * radius + center[0], center[1] - height, Math.sin(angle) * radius + center[2]));
		if (j > 0) {
			triangle(top[top.length-1], top[0], top[top.length-2]);
			triangle(bot[0], bot[bot.length-1], bot[bot.length-2]);
			corners.push(top[top.length-1]);
			corners.push(bot[bot.length-1]);
		}
	}
	for (let i = 0; i < num*2; i++) {
		quad((0+i)%(2*num), (1+i)%(2*num), (3+i)%(2*num), (2+i)%(2*num), corners);
	}
}

function MakePrism(A, B, C, D, height) {
	var corners = [
		A,
		B,
		C,
		D
	];
	let len = corners.length;
	for (let i = 0; i < len; i++) {
		corners.push(vec4(corners[i][0], corners[i][1] + height, corners[i][2], 1));
	}
	quad(0, 1, 2, 3, corners);
	quad(4, 5, 6, 7, corners);

	quad(0, 1, 5, 4, corners);
	quad(3, 2, 6, 7, corners);
	quad(1, 5, 6, 2, corners);
	quad(0, 3, 7, 4, corners);
}


function DrawLeg(root, which) {
	//The root represents one of the top corners of the legs and corresponds to an outer corner of the tabletop
	quad(0, 2, 6, 4); //AEGC
	quad(1, 3, 7, 5); //BDFH
	quad(0, 2, 3, 1); //ABCD
	quad(2, 3, 7, 6); //CDGH
	quad(4, 5, 7, 6); //EFGH
	quad(0, 1, 5, 4); //ABEF
}

var texCoord = [
    vec2(0, 0),
    vec2(0, 1),
    vec2(1, 1),
    vec2(1, 0)];

function quad(a, b, c, d, top) {
	var t1 = subtract(top[b], top[a]);
	var t2 = subtract(top[c], top[b]);
	var normal = normalize(vec3(cross(t1, t2)));

	pointsArray.push(top[a]);
	pointsArray.push(top[b]);
	pointsArray.push(top[c]);
	pointsArray.push(top[a]);
	pointsArray.push(top[c]);
	pointsArray.push(top[d]);
	for(let i = 0; i < 6; i++){
		normalsArray.push(normal);
	}

	texCoordsArray.push(texCoord[0]);
	texCoordsArray.push(texCoord[1]);
	texCoordsArray.push(texCoord[2]);
	texCoordsArray.push(texCoord[0]);
	texCoordsArray.push(texCoord[2]);
	texCoordsArray.push(texCoord[3]);

}

//---------------------------------------------------------------------------------------------------------------//


window.onload = function init()
{
    canvas = document.getElementById( "gl-canvas" );

    gl = WebGLUtils.setupWebGL( canvas );
    if ( !gl ) { alert( "WebGL isn't available" ); }

    gl.viewport( 0, 0, canvas.width, canvas.height );
    gl.clearColor( 1.0, 1.0, 1.0, 1.0 );

    gl.enable(gl.DEPTH_TEST);

    //
    //  Load shaders and initialize attribute buffers
    //
    program = initShaders( gl, "vertex-shader", "fragment-shader" );
    gl.useProgram( program );

	  openNewTexture('sky.jpg');
	  openNewTexture('floor.jpg');
	  openNewTexture('needles.jpg');
	  openNewTexture('gold.jpg');
	  openNewTexture('wall.jpg');
	  openNewTexture('wrapping.jpg');
    openNewTexture('rug.jpg');

  	MakePresent();
    MakeTable();
  	MakeStar();
  	MakeOrnament();
  	MakeGingerbreadMan();
  	MakeChristmasTree();
  	MakeMilk();
  	MakeCookies();
  	MakeRoom();

    var nBuffer = gl.createBuffer();
    gl.bindBuffer( gl.ARRAY_BUFFER, nBuffer );
    gl.bufferData( gl.ARRAY_BUFFER, flatten(normalsArray), gl.STATIC_DRAW );

    var vNormal = gl.getAttribLocation( program, "vNormal" );
    gl.vertexAttribPointer( vNormal, 3, gl.FLOAT, false, 0, 0 );
    gl.enableVertexAttribArray( vNormal );

    var vBuffer = gl.createBuffer();
    gl.bindBuffer( gl.ARRAY_BUFFER, vBuffer );
    gl.bufferData( gl.ARRAY_BUFFER, flatten(pointsArray), gl.STATIC_DRAW );

    var vPosition = gl.getAttribLocation(program, "vPosition");
    gl.vertexAttribPointer(vPosition, 4, gl.FLOAT, false, 0, 0);
    gl.enableVertexAttribArray(vPosition);

    thetaLoc = gl.getUniformLocation(program, "theta");

    projectionMatrix = ortho(-8, 8, -8, 8, -8, 8);

	modelViewLoc = gl.getUniformLocation( program, "modelViewMatrix" );
    projectionMatrixLoc = gl.getUniformLocation( program, "projectionMatrix" );

	gl.uniformMatrix4fv( gl.getUniformLocation(program, "projectionMatrix"),
       false, flatten(projectionMatrix));


	// These four just set the handlers for the buttons.
	document.getElementById("thetaup").addEventListener("click", function(e) {
		AllInfo.theta += AllInfo.dr;
	});
	document.getElementById("thetadown").addEventListener("click", function(e) {
		AllInfo.theta -= AllInfo.dr;
	});
	document.getElementById("phiup").addEventListener("click", function(e) {
		AllInfo.phi += AllInfo.dr;
	});
	document.getElementById("phidown").addEventListener("click", function(e) {
		AllInfo.phi -= AllInfo.dr;
	});

	// Set the scroll wheel to change the zoom factor.
	document.getElementById("gl-canvas").addEventListener("wheel", function(e) {
	 if (e.wheelDelta > 0) {
		 AllInfo.zoomFactor = Math.max(0.1, AllInfo.zoomFactor - 0.1);
	 } else {
		 AllInfo.zoomFactor += 0.1;
	 }
	});

	//************************************************************************************
	//* When you click a mouse button, set it so that only that button is seen as
	//* pressed in AllInfo. Then set the position. The idea behind this and the mousemove
	//* event handler's functionality is that each update we see how much the mouse moved
	//* and adjust the camera value by that amount.
	//************************************************************************************
	document.getElementById("gl-canvas").addEventListener("mousedown", function(e) {
	 if (e.which == 1) {
		 AllInfo.mouseDownLeft = true;
		 AllInfo.mouseDownRight = false;
		 AllInfo.mousePosOnClickY = e.y;
		 AllInfo.mousePosOnClickX = e.x;
	 } else if (e.which == 3) {
		 AllInfo.mouseDownRight = true;
		 AllInfo.mouseDownLeft = false;
		 AllInfo.mousePosOnClickY = e.y;
		 AllInfo.mousePosOnClickX = e.x;
	 }
	});

	document.addEventListener("mouseup", function(e) {
	 AllInfo.mouseDownLeft = false;
	 AllInfo.mouseDownRight = false;
	});

	document.addEventListener("mousemove", function(e) {
	 if (AllInfo.mouseDownRight) {
		 AllInfo.translateX += (e.x - AllInfo.mousePosOnClickX)/90;
		 AllInfo.mousePosOnClickX = e.x;

		 AllInfo.translateY -= (e.y - AllInfo.mousePosOnClickY)/90;
		 AllInfo.mousePosOnClickY = e.y;
	 } else if (AllInfo.mouseDownLeft) {
		 AllInfo.phi += (e.x - AllInfo.mousePosOnClickX)/100;
		 AllInfo.mousePosOnClickX = e.x;

		 AllInfo.theta += (e.y - AllInfo.mousePosOnClickY)/100;
		 AllInfo.mousePosOnClickY = e.y;
	 }
	});

	var tBuffer = gl.createBuffer();
    gl.bindBuffer( gl.ARRAY_BUFFER, tBuffer );
    gl.bufferData( gl.ARRAY_BUFFER, flatten(texCoordsArray), gl.STATIC_DRAW );

    var vTexCoord = gl.getAttribLocation( program, "vTexCoord" );
	console.log(vTexCoord);
    gl.vertexAttribPointer(vTexCoord, 2, gl.FLOAT, false, 0, 0 );
    gl.enableVertexAttribArray(vTexCoord );

    render();
}

document.onkeyup = checkKeyUp;

function checkKeyUp(e) {
	 e = e || window.event;
	if (e.keyCode == '65') {
		rotate_star = !rotate_star
	}
}


function openNewTexture(picName)
{
    var i = textures.length;
    textures[i] = gl.createTexture();
    textures[i].image = new Image();
    textures[i].image.src = picName;
    textures[i].image.onload = function() { loadNewTexture(i); }
}

function loadNewTexture(index)
{
    gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
    gl.activeTexture(gl.TEXTURE0 + index);
    gl.bindTexture(gl.TEXTURE_2D, textures[index]);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGB, gl.RGB, gl.UNSIGNED_BYTE, textures[index].image);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.MIRRORED_REPEAT);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.MIRRORED_REPEAT);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
}

function setUp(a, d, s, m) {
	gl.uniform4fv(gl.getUniformLocation(program, "lightPosition"), flatten(lightPosition));

	gl.uniform4fv(gl.getUniformLocation(program, "ambientProduct"), flatten(a));
    gl.uniform4fv(gl.getUniformLocation(program, "diffuseProduct"), flatten(d));
    gl.uniform4fv(gl.getUniformLocation(program, "specularProduct"), flatten(s));
    gl.uniform1f(gl.getUniformLocation(program, "shininess"), m);
}

function drawTable(x, y, z) {
	materialAmbient = vec4( .35, .2, .2, 1.0 );
	materialDiffuse = vec4( .8, .4, .4, 1.0);
	materialSpecular = vec4( .25, .25, .25, 1 );
	materialShininess = 50.0;
	let ambientProduct = mult(lightAmbient, materialAmbient);
    let diffuseProduct = mult(lightDiffuse, materialDiffuse);
    let specularProduct = mult(lightSpecular, materialSpecular);
	setUp(ambientProduct, diffuseProduct, specularProduct, materialShininess);

	modelView = lookAt(eye, at, up);

	modelView = mult(modelView, translate(x, y, z));
	modelView = mult(modelView, scale4(1.05, 1.05, 1.05));

    gl.uniformMatrix4fv( gl.getUniformLocation(program,
            "modelViewMatrix"), false, flatten(modelView));

	gl.drawArrays(gl.TRIANGLES, tableIndex, tableNum);
}

function drawPlate(x, y, z) {
	//drawplate
	let materialAmbient = vec4( .4, .4, .6, 1 );
	let materialDiffuse = vec4( 0.4, .4, 0.3, 1);
	let materialSpecular = vec4( .2, .2, .1, 1 );
	let materialShininess = 25.0;

	let ambientProduct = mult(lightAmbient, materialAmbient);
	let diffuseProduct = mult(lightDiffuse, materialDiffuse);
	let specularProduct = mult(lightSpecular, materialSpecular);
	setUp(ambientProduct, diffuseProduct, specularProduct, materialShininess);

	modelView = lookAt(eye, at, up);

	modelView = mult(modelView, translate(x, y, z));
	modelView = mult(modelView, scale4(.15, .15, .15));
	gl.uniformMatrix4fv( gl.getUniformLocation(program,
					"modelViewMatrix"), false, flatten(modelView) );

	gl.drawArrays( gl.TRIANGLES, plateIndex, plateNum);
}

function drawCookies(x, y, z) {

	//draw cookies
	materialAmbient = vec4( .4, .3, .2, 1.0 );
	materialDiffuse = vec4( .4, .35, .2, 1.0);
	materialSpecular = vec4( .4, .35, .2, 1 );
	materialShininess = 100.0;

	ambientProduct = mult(lightAmbient, materialAmbient);
	diffuseProduct = mult(lightDiffuse, materialDiffuse);
	specularProduct = mult(lightSpecular, materialSpecular);
	setUp(ambientProduct, diffuseProduct, specularProduct, materialShininess);

	modelView = lookAt(eye, at, up);

	modelView = mult(modelView, translate(x, y, z));
	modelView = mult(modelView, scale4(.05, .05, .05));
	gl.uniformMatrix4fv( gl.getUniformLocation(program,
					"modelViewMatrix"), false, flatten(modelView) );

	gl.drawArrays( gl.TRIANGLES, cookieIndex, cookieNum);

	materialAmbient = vec4( .05, .05, 0, 1 );
	materialDiffuse = vec4( 0.2, .2, 0., 1);
	materialSpecular = vec4( .2, .2, .0, 1 );
	materialShininess = 25.0;

	ambientProduct = mult(lightAmbient, materialAmbient);
	diffuseProduct = mult(lightDiffuse, materialDiffuse);
	specularProduct = mult(lightSpecular, materialSpecular);
	setUp(ambientProduct, diffuseProduct, specularProduct, materialShininess);

	modelView = lookAt(eye, at, up);

	modelView = mult(modelView, translate(x, y, z));
	modelView = mult(modelView, scale4(.05, .05, .05));
	gl.uniformMatrix4fv( gl.getUniformLocation(program,
					"modelViewMatrix"), false, flatten(modelView) );

	gl.drawArrays( gl.TRIANGLES, chipIndex, chipNum);
}

function drawPresent(x, y, z, s, w) {
	materialAmbient = vec4( 1, 1, 1, 1.0 );
	materialDiffuse = vec4( 1, 1, 1, 1.0);
	materialSpecular = vec4( .25, .25, .25, 1 );
	materialShininess = 50.0;
	let ambientProduct = mult(lightAmbient, materialAmbient);
		let diffuseProduct = mult(lightDiffuse, materialDiffuse);
		let specularProduct = mult(lightSpecular, materialSpecular);
	setUp(ambientProduct, diffuseProduct, specularProduct, materialShininess);

	modelView = lookAt(eye, at, up);

	modelView = mult(modelView, translate(x, y, z));
	modelView = mult(modelView, scale4(s, w, s));
	gl.uniform1i(gl.getUniformLocation(program, "textureFlag"), true);
		gl.uniformMatrix4fv( gl.getUniformLocation(program,
						"modelViewMatrix"), false, flatten(modelView));
	  gl.drawArrays( gl.TRIANGLES, presentIndex, presentNum);
}

function drawStar(x, y, z, rot) {
	materialAmbient = vec4( .85, .85, 0., 1.0 );
	materialDiffuse = vec4( .6, .6, 0., 1.0);
	materialSpecular = vec4( .8, .8, 0, 1 );
	materialShininess = 50.0;
	let ambientProduct = mult(lightAmbient, materialAmbient);
    let diffuseProduct = mult(lightDiffuse, materialDiffuse);
    let specularProduct = mult(lightSpecular, materialSpecular);
	setUp(ambientProduct, diffuseProduct, specularProduct, materialShininess);

	modelView = lookAt(eye, at, up);

    modelView = mult(modelView, rotate(rot, [0, 1, 0] ));

	modelView = mult(modelView, translate(x, y, z));

    gl.uniformMatrix4fv( gl.getUniformLocation(program,
            "modelViewMatrix"), false, flatten(modelView));

	gl.drawArrays(gl.TRIANGLES, starIndex, starNum);
}

function drawGingerbreadMan(x, y, z) {
	materialAmbient = vec4( .35, .2, .15, 1.0 );
	materialDiffuse = vec4( .4, .2, .2, 1.0);
	materialSpecular = vec4( .5, .25, .25, 1 );
	materialShininess = 25.0;

	let ambientProduct = mult(lightAmbient, materialAmbient);
	let diffuseProduct = mult(lightDiffuse, materialDiffuse);
	let specularProduct = mult(lightSpecular, materialSpecular);
	setUp(ambientProduct, diffuseProduct, specularProduct, materialShininess);

	modelView = lookAt(eye, at, up);

	modelView = mult(modelView, translate(x, y, z));
	modelView = mult(modelView, scale4(.2, .2, .2));
	gl.uniformMatrix4fv(modelViewLoc, false, flatten(modelView) );
	gl.uniformMatrix4fv(projectionMatrixLoc, false, flatten(projectionMatrix) );
	gl.uniformMatrix4fv( gl.getUniformLocation(program,
		"modelViewMatrix"), false, flatten(modelView));

	gl.drawArrays(gl.TRIANGLES, gingerbreadIndex, gingerbreadNum);
}

function drawOrnament(ox, oy, oz) {
	materialAmbient = vec4( .6, .8, .8, 1.0 );
	materialDiffuse = vec4( .5, .5, .5, 1.0);
	materialSpecular = vec4( .5, .5, .5, 1 );
	materialShininess = 10.0;

	let ambientProduct = mult(lightAmbient, materialAmbient);
	let diffuseProduct = mult(lightDiffuse, materialDiffuse);
	let specularProduct = mult(lightSpecular, materialSpecular);
	setUp(ambientProduct, diffuseProduct, specularProduct, materialShininess);

	modelView = lookAt(eye, at, up);

	modelView = mult(modelView, translate(-ox, -oy, oz));
	modelView = mult(modelView, scale4(.1, .1, .1));
	//First Sphere
	modelView = mult(modelView, translate(0,-.5,1,0));
	gl.uniformMatrix4fv(modelViewLoc, false, flatten(modelView) );
	gl.uniformMatrix4fv(projectionMatrixLoc, false, flatten(projectionMatrix) );
	gl.uniformMatrix4fv( gl.getUniformLocation(program,
		"modelViewMatrix"), false, flatten(modelView));
	gl.drawArrays(gl.TRIANGLES, ornamentIndex, ornamentNum);

		//Second Sphere
	modelView = mult(modelView, translate(0,1.4,0,0));
	gl.uniformMatrix4fv(modelViewLoc, false, flatten(modelView) );
	gl.uniformMatrix4fv(projectionMatrixLoc, false, flatten(projectionMatrix) );
	gl.drawArrays(gl.TRIANGLES, ornamentIndex, ornamentNum);

	gl.uniformMatrix4fv(modelViewLoc, false, flatten(modelView) );
	gl.uniformMatrix4fv(projectionMatrixLoc, false, flatten(projectionMatrix) );
	gl.drawArrays(gl.TRIANGLES, topIndex, topNum);

}
function drawChristmasTree(x, y, z)
{
	var materialAmbient = vec4( 0.1, .3, .1, 1.0 );
	var materialDiffuse = vec4( 0.5, .5, 0.5, 1.0);
	var materialSpecular = vec4( .3, .5, .3, 1.0 );
	var materialShininess = 10.0;

	let ambientProduct = mult(lightAmbient, materialAmbient);
	let diffuseProduct = mult(lightDiffuse, materialDiffuse);
	let specularProduct = mult(lightSpecular, materialSpecular);
	setUp(ambientProduct, diffuseProduct, specularProduct, materialShininess);

	modelView = lookAt(eye, at, up);


	modelView = mult(modelView, translate(x, y, z));
	modelView = mult(modelView, scale4(5, -5, 5));
	gl.uniformMatrix4fv( gl.getUniformLocation(program,
					"modelViewMatrix"), false, flatten(modelView) );

	//gl.drawArrays( gl.TRIANGLES, 0, 24*6);
	gl.drawArrays( gl.TRIANGLES, treeIndex, treeNum);
}


function drawRoom(x, y, z) {
	gl.uniform1i(gl.getUniformLocation(program, "textureFlag"), true);
	gl.uniform1i(gl.getUniformLocation(program, "texture"), 1);
	materialAmbient = vec4( .4, .3, .3, 1.0 );
	materialDiffuse = vec4( .5, .5, .5, 1.0);
	materialSpecular = vec4( .3, .3, .3, 1 );
	materialShininess = 50.0;
	var materialShininess = 100.0;

	let ambientProduct = mult(lightAmbient, materialAmbient);
	let diffuseProduct = mult(lightDiffuse, materialDiffuse);
	let specularProduct = mult(lightSpecular, materialSpecular);
	setUp(ambientProduct, diffuseProduct, specularProduct, materialShininess);

	modelView = lookAt(eye, at, up);

	modelView = mult(modelView, translate(-.75, -.1, -.15));

	modelView = mult(modelView, translate(x, y, z));
	gl.uniformMatrix4fv( gl.getUniformLocation(program,
					"modelViewMatrix"), false, flatten(modelView) );

	gl.drawArrays( gl.TRIANGLES, floorIndex, floorNum);

	gl.uniform1i(gl.getUniformLocation(program, "texture"), 4);
	materialAmbient = vec4( .55, .55, .65, 1 );
	materialDiffuse = vec4( 0.4, .4, 0.3, 1);
	materialSpecular = vec4( .2, .2, .1, 1 );
	materialShininess = 100.0;

	ambientProduct = mult(lightAmbient, materialAmbient);
	diffuseProduct = mult(lightDiffuse, materialDiffuse);
	specularProduct = mult(lightSpecular, materialSpecular);
	setUp(ambientProduct, diffuseProduct, specularProduct, materialShininess);
	gl.drawArrays( gl.TRIANGLES, wallIndex, wallNum);
}

function drawRug(x,y,z)
{
  materialAmbient = vec4( .8, .8, .8, 1.0 );
  materialDiffuse = vec4( .6, .5, .5, 1.0);
  materialSpecular = vec4( .3, .3, .3, 1 );
  materialShininess = 50.0;
  var materialShininess = 100.0;

  let ambientProduct = mult(lightAmbient, materialAmbient);
  let diffuseProduct = mult(lightDiffuse, materialDiffuse);
  let specularProduct = mult(lightSpecular, materialSpecular);
  setUp(ambientProduct, diffuseProduct, specularProduct, materialShininess);

  modelView = lookAt(eye, at, up);

  modelView = mult(modelView, translate(-.75, -.7, -.15));

  modelView = mult(modelView, translate(x, y, z));
  	modelView = mult(modelView, scale4(.7, .7, .7));

  gl.uniformMatrix4fv( gl.getUniformLocation(program,
          "modelViewMatrix"), false, flatten(modelView) );

  gl.drawArrays( gl.TRIANGLES, floorIndex, floorNum);
}

function drawMilk(x, y, z)
{
	var materialAmbient = vec4( .4, .4, .6, 1 );
	var materialDiffuse = vec4( 0.4, .4, 0.3, 1);
	var materialSpecular = vec4( .2, .2, .1, 1 );
	var materialShininess = 100.0;

	let ambientProduct = mult(lightAmbient, materialAmbient);
	let diffuseProduct = mult(lightDiffuse, materialDiffuse);
	let specularProduct = mult(lightSpecular, materialSpecular);
	setUp(ambientProduct, diffuseProduct, specularProduct, materialShininess);

	modelView = lookAt(eye, at, up);

	modelView = mult(modelView, translate(-.75, -.1, -.15));

	modelView = mult(modelView, translate(x, y, z));
	gl.uniformMatrix4fv( gl.getUniformLocation(program,
					"modelViewMatrix"), false, flatten(modelView) );

	//gl.drawArrays( gl.TRIANGLES, 0, 24*6);
	gl.drawArrays( gl.TRIANGLES, cupIndex, cupNum);

	materialAmbient = vec4( .95, .95, .85, 1 );
	materialDiffuse = vec4( 0.4, .4, 0.3, 1);
	materialSpecular = vec4( .2, .2, .1, 1 );
	materialShininess = 100.0;

	ambientProduct = mult(lightAmbient, materialAmbient);
	diffuseProduct = mult(lightDiffuse, materialDiffuse);
	specularProduct = mult(lightSpecular, materialSpecular);
	setUp(ambientProduct, diffuseProduct, specularProduct, materialShininess);
	gl.drawArrays( gl.TRIANGLES, milkIndex, milkNum);

}

document.onkeyup = checkKeyUp;
var audio = document.createElement('audio');
audio.src = 'winterwonder.mp3';
audio.currentTime = 0;

function checkKeyUp(e) {
	 e = e || window.event;
	if (e.keyCode == '65') {
		rotate_star = !rotate_star
		if(rotate_star)
		{
			audio.play();
		}
		else {
			audio.pause();
		}
	}
}

var a = .25;
var b = .9;

var c = -.25;
var d = 1.3;
var w = -.6;
//done
var e = -.8;
var f = 1.9;
//done
var g = .65;
var h = 2.3;

var rotation = 0;

// Variables that control the orthographic projection bounds.
var y_max = 3;
var y_min = -3;
var x_max = 4;
var x_min = -4;
var near = -10;
var far = 10;

function render() {
    gl.clear( gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

	projectionMatrix = ortho( x_min*AllInfo.zoomFactor - AllInfo.translateX,
														x_max*AllInfo.zoomFactor - AllInfo.translateX,
														y_min*AllInfo.zoomFactor - AllInfo.translateY,
														y_max*AllInfo.zoomFactor - AllInfo.translateY,
														near, far);
	gl.uniformMatrix4fv(projectionMatrixLoc, false, flatten(projectionMatrix));

	eye = vec3( AllInfo.radius*Math.cos(AllInfo.phi),
							AllInfo.radius*Math.sin(AllInfo.theta),
							AllInfo.radius*Math.sin(AllInfo.phi));

	if (rotate_star) {
		rotation = rotation += 1.0;
	}

	gl.uniform1i(gl.getUniformLocation(program, "textureFlag"), true);
	gl.uniform1i(gl.getUniformLocation(program, "texture"), 1);
	drawTable(6, -.45, -.5);
	drawGingerbreadMan(3, -1.31, .3);
	gl.uniform1i(gl.getUniformLocation(program, "texture"), 2);
	drawChristmasTree(0, .2, 0);
	gl.uniform1i(gl.getUniformLocation(program, "textureFlag"), false);
	drawOrnament(a, b, a);
	drawOrnament(c, d, w);
	drawOrnament(e, f, c);
	drawOrnament(g, h, g);
	gl.uniform1i(gl.getUniformLocation(program, "textureFlag"), true);
	gl.uniform1i(gl.getUniformLocation(program, "texture"), 3);
	drawStar(3, 0, 0, rotation);
	gl.uniform1i(gl.getUniformLocation(program, "textureFlag"), false);
	drawMilk(6, -.5, -.6);
	gl.uniform1i(gl.getUniformLocation(program, "textureFlag"), true);
	gl.uniform1i(gl.getUniformLocation(program, "texture"), 0);
	drawPlate(3, -1.31, .1);
	gl.uniform1i(gl.getUniformLocation(program, "textureFlag"), false);
	drawCookies(5-2.1, -1.2, -.1);
	drawCookies(5-1.8, -1.2, .1);
	drawCookies(5-2.2, -1.2, .2);
	gl.uniform1i(gl.getUniformLocation(program, "textureFlag"), true);
	gl.uniform1i(gl.getUniformLocation(program, "texture"), 5);
	drawPresent(1, -2.1, 1, .7, .7);
	drawPresent(0, -2.1, 1.5, .7, .5);
	drawPresent(-1, -2.1, 1.5, .5, .7);
  gl.uniform1i(gl.getUniformLocation(program, "textureFlag"), true);
	gl.uniform1i(gl.getUniformLocation(program, "texture"), 6);
  drawRug(3, -1, 2);
	drawRoom(3.3, -1.4, 2.9);

    requestAnimFrame(render);
}
