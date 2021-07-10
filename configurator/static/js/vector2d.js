/*
Simple 2D JavaScript Vector2d Class
Hacked from evanw's lightgl.js
https://github.com/evanw/lightgl.js/blob/master/src/vector.js
*/

function Vector2d(x, y) {
	this.x = x || 0;
	this.y = y || 0;
}

/* INSTANCE METHODS */

Vector2d.prototype = {
	negative: function() {
		this.x = -this.x;
		this.y = -this.y;
		return this;
	},
	add: function(v) {
		if (v instanceof Vector2d) {
			this.x += v.x;
			this.y += v.y;
		} else {
			this.x += v;
			this.y += v;
		}
		return this;
	},
	subtract: function(v) {
		if (v instanceof Vector2d) {
			this.x -= v.x;
			this.y -= v.y;
		} else {
			this.x -= v;
			this.y -= v;
		}
		return this;
	},
	multiply: function(v) {
		if (v instanceof Vector2d) {
			this.x *= v.x;
			this.y *= v.y;
		} else {
			this.x *= v;
			this.y *= v;
		}
		return this;
	},
	divide: function(v) {
		if (v instanceof Vector2d) {
			if(v.x != 0) this.x /= v.x;
			if(v.y != 0) this.y /= v.y;
		} else {
			if(v != 0) {
				this.x /= v;
				this.y /= v;
			}
		}
		return this;
	},
	equals: function(v) {
		return this.x == v.x && this.y == v.y;
	},
	dot: function(v) {
		return this.x * v.x + this.y * v.y;
	},
	cross: function(v) {
		return this.x * v.y - this.y * v.x
	},
	length: function() {
		return Math.sqrt(this.dot(this));
	},
	normalize: function() {
		return this.divide(this.length());
	},
	min: function() {
		return Math.min(this.x, this.y);
	},
	max: function() {
		return Math.max(this.x, this.y);
	},
	toAngles: function() {
		return -Math.atan2(-this.y, this.x);
	},
	angleTo: function(a) {
		return Math.acos(this.dot(a) / (this.length() * a.length()));
	},
	toArray: function(n) {
		return [this.x, this.y].slice(0, n || 2);
	},
	clone: function() {
		return new Vector2d(this.x, this.y);
	},
	set: function(x, y) {
		this.x = x; this.y = y;
		return this;
	}
};

/* STATIC METHODS */
Vector2d.fromAngles = function(alpha) {
  return new Vector2d(Math.cos(alpha), Math.sin(alpha));
};
Vector2d.randomDirection = function() {
  return Vector2d.fromAngles(Math.random() * Math.PI * 2);
};
Vector2d.negative = function(v) {
	return new Vector2d(-v.x, -v.y);
};
Vector2d.add = function(a, b) {
	if (b instanceof Vector2d) return new Vector2d(a.x + b.x, a.y + b.y);
	else return new Vector2d(a.x + b, a.y + b);
};
Vector2d.subtract = function(a, b) {
	if (b instanceof Vector2d) return new Vector2d(a.x - b.x, a.y - b.y);
	else return new Vector2d(a.x - b, a.y - b);
};
Vector2d.multiply = function(a, b) {
	if (b instanceof Vector2d) return new Vector2d(a.x * b.x, a.y * b.y);
	else return new Vector2d(a.x * b, a.y * b);
};
Vector2d.divide = function(a, b) {
	if (b instanceof Vector2d) return new Vector2d(a.x / b.x, a.y / b.y);
	else return new Vector2d(a.x / b, a.y / b);
};
Vector2d.equals = function(a, b) {
	return a.x == b.x && a.y == b.y;
};
Vector2d.dot = function(a, b) {
	return a.x * b.x + a.y * b.y;
};
Vector2d.cross = function(a, b) {
	return a.x * b.y - a.y * b.x;
};