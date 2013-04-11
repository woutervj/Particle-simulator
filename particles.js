var N = 343;
var U;
var d = 20;
var dw = 10;
var d2 = d*d;
var fieldWidth = 800;
var fieldHeight = 800;
var G= 10;
var k=G/d;
var v0 = 0;
var g = 40;
var core = Math.pow(2,1/6)*d;

var vmax = fieldWidth/10;
var Fmax = 50;

function particle(m,x,y,vx,vy)
{
	this.m=m;
	this.x=x;
	this.y=y;
	this.vx=vx;
	this.vy=vy;
	this.Fx = 0;
	this.Fy = 0;
	this.color = "red";
	this.d = d;

	this.step1 = function (dt,u) {
		F = u.F(this);
		var ax = F.x/this.m;
		var ay = F.y/this.m;
		this.vx += ax * dt; this.vx = Math.max(-vmax, Math.min(vmax,this.vx));
		this.vy += ay * dt; this.vy = Math.max(-vmax, Math.min(vmax,this.vy));
		//this.vx *=0.5;
		//this.vy *=0.5;
		this.x += this.vx*dt;
		this.y += this.vy*dt;
		this.checkBounds(u.w, u.h);
	};
	
	this.step = function (dt,u) {
		var ax = this.Fx/this.m;
		var ay = this.Fy/this.m;
		this.vx += ax * dt; this.vx = Math.max(-vmax, Math.min(vmax,this.vx));
		this.vy += ay * dt; this.vy = Math.max(-vmax, Math.min(vmax,this.vy));
		this.x += this.vx*dt;
		this.y += this.vy*dt;
		this.checkBounds(u.w, u.h);
	};
	
	
	this.draw = function(ctx)
	{
		ctx.beginPath();
		ctx.arc(this.x, this.y, dw, 0, 2*Math.PI);
		ctx.fillStyle=this.color;
		ctx.fill();
	};
	
	this.checkBounds = function (w,h) {
		if (this.x <dw) {
			this.x = 2*dw - this.x;
			this.vx = - this.vx;
		}
		if (this.y <dw) {
			this.y = 2*dw- this.y;
			this.vy = - this.vy;
		}
		if (this.x>w) {
			this.x = 2*w - this.x;
			this.vx = -this.vx;
		}
		if (this.y>h) {
			this.y = 2*h - this.y;
			this.vy = -this.vy;
		}
	};
	
	this.checkBoundsWrap = function (w,h) {
		if (this.x < 0) {
			this.x+=w;
		}
		if (this.y < 0) {
			this.y += h;
		}
		if (this.x>w) {
			this.x -= w;
		}
		if (this.y>h) {
			this.y -= h;
		}
	};
	
	this.clearForce = function () {
		this.Fx = 0;
		this.Fy = 0;
	};
}

function universe(w,h)
{
	this.particles = new Array;
	this.addParticle = function (aParticle) { this.particles.push(aParticle); };
	this.dt = 0.05;
	this.running = false;
	this.w = w;
	this.h = h;
	this.canvas=null;
	this.fields=new Array();
	this.drawfields=new Array();
	
	this.addField = function (field) {
		//add the function to compute the field function (point, universe)
		this.fields.push(field.func);
		//add the variables associated with the field to allow computation and drawing
		//should be prefixed with the field name.
		for (var key in field.vars) {		
			this[key] = field.vars[key];
		}
		//add a function to draw the field. function (context, universe)
		this.drawfields.push(field.draw);
	}

	this.F = function (aParticle) {
		var Fx = 0;
		var Fy = 0;
		var F;
		for (var i=0; i<this.particles.length; i++) {
			if (this.particles[i] != aParticle) {
				F = this.force(aParticle,this.particles[i]); 
				Fx += F.x;
				Fy += F.y;
			}
		}
		F = this.fieldForce(aParticle);
		Fx += F.x;
		Fy += F.y;		
		return {x: Fx, y: Fy};
	};
	
	this.force = LJforce; 
	
	this.fieldForce = function (p) {
		var Fx=0, Fy=0;
		for (var i=0; i< this.fields.length; i++) {
			var m = this; 
			var F = this.fields[i](p,m);
			Fx += F.x;
			Fy += F.y;
		}
		return {x:Fx, y:Fy};
	};
	
	this.step1 = function () {
		for (var i=0; i< this.particles.length; i++) {
			this.particles[i].step(this.dt,this);
			this.particles[i].update();
		}
	};
	
	this.step = function () {
		//clear all the forces
		for (var i=0; i< this.particles.length; i++) {
			this.particles[i].clearForce();
		}
		//compute the new interaction forces
		for (var i=0; i< this.particles.length; i++) {
			for (var j=i+1; j< this.particles.length; j++) {
				var F = this.force(this.particles[i],this.particles[j]);
				// apply Newton's third law
				this.particles[i].Fx += F.x;
				this.particles[i].Fy += F.y;
				this.particles[j].Fx -= F.x;
				this.particles[j].Fy -= F.y;
			}
		}
		
		//compute force fields
		for (var i=0; i< this.particles.length; i++) {
			var F = this.fieldForce(this.particles[i]);
				this.particles[i].Fx += F.x;
				this.particles[i].Fy += F.y;			
		}
		
		//and let the particles go!
		for (var i=0; i< this.particles.length; i++) {
			this.particles[i].step(this.dt,this);
		}
		this.draw();		
	};
	
	this.draw = function ()
	{
		var ctx = this.canvas.getContext("2d");
		ctx.clearRect(0,0,fieldWidth,fieldHeight);
		for (var i=0; i<this.drawfields.length; i++) {
			var m = this;
			this.drawfields[i](ctx, m);
		}
		for (var i=0; i< this.particles.length; i++) {
			this.particles[i].draw(ctx);
		}
	}
	
	this.go = function()
	{
		this.running = true;
		this.timedStep()
	};
	
	this.timedStep = function () {
        var mdl = this;
        window.setTimeout(function () { mdl.step(); if (mdl.running) {mdl.timedStep();} }, 10); 
    };
    
    this.stop = function ()
    {
    	this.running = false;
    }
}

function init()
{
	U = new universe(fieldWidth, fieldHeight);
	U.canvas = document.getElementById('playground');
	U.canvas.width = U.w;
	U.canvas.height = U.h;
	U.addField({func:LJBox, draw: function (c,u) {}, vars: {}});
	U.addField({func:funnelForce, draw: funnelDraw, vars: {funnelLines: funnel}});
	var px = 2*d; var py = d*1.5;
	for (var i=0; i<N; i++) {
		var fw = fieldWidth - 2*dw;
		var p = new particle(1,px,py,2*v0*Math.random()-v0,2*v0*Math.random()-v0);
		px += d;
		if (px>fw-py) {py +=d; px = py+0.5*d; }
		U.addParticle(p);
		U.draw();
	}
}

//these functions come in useful
function sqr(x) { return x * x }
function dist2(v, w) { return sqr(v.x - w.x) + sqr(v.y - w.y); }
function length2(p) {return sqr(p.x) + sqr(p.y);}
function length(p) {return Math.sqrt(length2(p));}
function normalize(p) {var l=length(p); return {x: p.x/l, y:p.y/l};}

function distToSegmentSquared(p, v, w) {
  var l2 = dist2(v, w);
  if (l2 == 0) return dist2(p, v);
  var t = ((p.x - v.x) * (w.x - v.x) + (p.y - v.y) * (w.y - v.y)) / l2;
  if (t < 0) return 10*d;//dist2(p, v);
  if (t > 1) return 10*d;//dist2(p, w);
  return dist2(p, { x: v.x + t * (w.x - v.x),
                    y: v.y + t * (w.y - v.y) });
}


function distToSegment(p, v, w) { return Math.sqrt(distToSegmentSquared(p, v, w));}


//functions defining the force fields
//Lennard-Jones with repulsive boundaries

function LJBox (p,u) {

	var Fx = 0;
	var Fy = p.m*g;

	if (p.x < 3*d) { 
		var delta = d/(2*p.x); 
		d6=delta*delta*delta; 
		d12 = d6*d6; 
		Fx += G*(12*d12 - 6*d6)/(2*p.x);
	};
	if (p.y < 3*d) { 
		var delta = d/p.y; 
		d6=delta*delta*delta; 
		d12 = d6*d6; 
		Fy += G*(12*d12 - 6*d6)/(2*p.x);
	};
	if (p.x > u.w-3*d) { 
		var delta = d/(2*(u.w-p.x)); 
		d6=delta*delta*delta; 
		d12 = d6*d6; 
		Fx += -G*(12*d12 - 6*d6)/(2*(u.w-p.x));
	};
	if (p.y > u.h-3*d) { 
		var delta = d/(2*(u.h-p.y)); 
		d6=delta*delta*delta; 
		d12 = d6*d6; 
		Fy += -G*(12*d12 - 6*d6)/(2*(u.h-p.y));
	};
	return {x: Fx, y: Fy};		
}

function sPath(points)
{
	this.points = new Array();
	this.lengthSquared = 0;
	for (var i=0; i<points.length-1; i+=2) {
		this.points.push({x:points[i], y:points[i+1]});
		if (i>0) this.lengthSquared += dist2(this.points[i/2], this.points[i/2-1])
	}	
			
	this.closestPointOnLine = function(p,i) {
		var p0 = this.points[i];
		var p1 = this.points[i+1];
		var l2 = dist2(p0,p1)
	  	if (l2 ==0) return this.p1;
	  	var t = ((p.x - p0.x) * (p1.x - p0.x) + (p.y - p0.y) * (p1.y - p0.y)) / l2;
	  	if (t < 0) return p0;	
	  	if (t > 1) return p1;
	  	return { x: p0.x + t * (p1.x - p0.x), y: p0.y + t * (p1.y - p0.y) };
	}
	
	this.closestPointOnPath = function (p)
	{
		var pc;
		var d = sqr(fieldWidth);
		for (var i=0; i<this.points.length-1; i++) {
			var pn = this.closestPointOnLine(p,i);
			var d2 = dist2(p,pn);
			if (d2 < d) {
				d = d2;
				pc = pn;
			}
		}
		return pc;
	}

	this.closestPointOnLine = function(p,i) {
		var p0 = this.points[i];
		var p1 = this.points[i+1];
		var l2 = dist2(p0,p1)
	  	if (l2 ==0) return this.p1;
	  	var t = ((p.x - p0.x) * (p1.x - p0.x) + (p.y - p0.y) * (p1.y - p0.y)) / l2;
	  	if (t < 0) return p0;	
	  	if (t > 1) return p1;
	  	return { x: p0.x + t * (p1.x - p0.x), y: p0.y + t * (p1.y - p0.y) };
	}

	
	this.parallelUnitVector = function (p)
	{
		var pc, uvx, uvy;
		var d = sqr(fieldWidth);
		var p0= this.points[0];
		var pc=p0;
		var pn=p0;
		var t;
		for (var i=1; i<this.points.length; i++) {
			var p1 = this.points[i];
			var l2 = dist2(p0,p1)
		  	if (l2 > 0) {
				t = ((p.x - p0.x) * (p1.x - p0.x) + (p.y - p0.y) * (p1.y - p0.y)) / l2;
				if (t < 0) pn = p0;	else if (t > 1) pn = p1;
				else {
					pn.x = p0.x + t * (p1.x - p0.x);
					pn.y = p0.y + t * (p1.y - p0.y);
				}		  	
		  	}
			var d2 = dist2(p,pn);
			if (d2 < d) {
				d = d2;
				pc = pn;
				uvx = p1.x-p0.x;
				uvy = p1.y-p0.y;
			}
		}
		return normalize({base: pc, direction: {x: uvx, y: uvy}}); 
	}
	
	this.draw = function(ctx)	{
		if (this.points.length <2) return;
		ctx.strokeStyle="grey";
		ctx.beginPath();
		ctx.moveTo(this.points[0].x,this.points[0].y);
		for (var i=1; i< this.points.length; i++) {
				ctx.lineTo(this.points[i].x,this.points[i].y);
		}
		ctx.stroke()			
	}
}


var fp = 1;
var funnel = [
	new sPath([0,0,
		fieldWidth/2-fp*d,fieldWidth/2-fp*d,
		fieldWidth/2-fp*d, fieldWidth/2+fp*d,
		0, fieldWidth/2+fp*d]),
	new sPath([fieldWidth,0,
		fieldWidth/2+fp*d, fieldWidth/2-fp*d,
		fieldWidth/2+fp*d, fieldWidth/2+fp*d,
		fieldWidth, fieldWidth/2+fp*d]),
	new sPath([2*fp*d, fieldWidth/2+3*fp*d,
		fieldWidth - 2*fp*d, fieldWidth/2+3*fp*d,
		fieldWidth - 2*fp*d, fieldWidth/2+4*fp*d,
		2*fp*d, fieldWidth/2+4*fp*d,
		2*fp*d, fieldWidth/2+3*fp*d])
	];

function funnelForce(p,u)
{
	var Fx=0, Fy=0;
	for (var i=0; i<u.funnelLines.length; i++) {
		var F = LJforce(p, u.funnelLines[i].closestPointOnPath(p));
		Fx += F.x;
		Fy += F.y;
	}
	return {x:Fx, y:Fy};
}

function funnelDraw(context,u)
{
	for (var i=0; i<u.funnelLines.length; i++) {
		u.funnelLines[i].draw(context);
	}
}

function LJforce (p1, p2) {
		var y = p2.y-p1.y;
		var x = p2.x-p1.x;
		var r2 = y*y + x*x;
		if (r2>7*d2) return {x:0,y:0};
		var r = Math.sqrt(r2);
		var delta = d2/r2;
		var d6 = delta*delta*delta;
		var d12 = d6*d6;
		var phi = Math.atan2(y,x);

		var F=-G*(12*d12 - 6*d6)/r;
		return {x: F*Math.cos(phi), y: F*Math.sin(phi)};
	};


function dragLineForce(p)
{
	
}

function funnelDraw(context,u)
{
	for (var i=0; i<u.draglines.length; i++) {
		u.dragLines[i].draw(context);
	}
}
