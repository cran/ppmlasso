decimalCount = function(x, max.dec = max(10, max(nchar(x))), tol = 1.e-1)
{
	digits  = 0:max.dec
	x.diff = (x - round(x, digits))/(10^(-1*(digits + 3)))
	num.dec = digits[min(which(abs(x.diff) < tol))]
	num.dec
}

spatRes = function(env.grid, coord = c("X", "Y"))
{
	x.col   = which(names(env.grid) == coord[1])
	y.col   = which(names(env.grid) == coord[2])
	x.uq    = sort(unique(env.grid[, x.col]))
	y.uq    = sort(unique(env.grid[, y.col]))
	n.dec   = max(unlist(lapply(x.uq, decimalCount)))
	x.diff  = diff(x.uq)
	y.diff  = diff(y.uq)
	x.dec   = unlist(lapply(x.diff, decimalCount))
	y.dec   = unlist(lapply(y.diff, decimalCount))
	x.step  = min(floor(x.diff*10^max(x.dec) + 0.1))/(10^max(x.dec))
	y.step  = min(floor(y.diff*10^max(y.dec) + 0.1))/(10^max(y.dec))
	return(list(x.step = x.step, y.step = y.step))
}

zapCoord = function(x, numdec = decimalCount(x))
{
	x.out  = floor(x*10^numdec + 0.1)/(10^numdec)
	x.out
}

griddify = function(envframe, tol = 0.01, coord = c("X", "Y"))
{
  xcol = match(coord[1], names(envframe))
  ycol = match(coord[2], names(envframe))
  quadx = envframe[,xcol]
  quady = envframe[,ycol]
  
  dx = quadx[duplicated(quadx)]
  dy = quady[duplicated(quady)]
  
  udx = sort(unique(dx))
  udy = sort(unique(dy))
  
  diffx = diff(udx)
  diffy = diff(udy)
  
  x_grain = median(diffx)
  y_grain = median(diffy)
  
  x_first = udx[which(diffx == x_grain)[1]]
  y_first = udy[which(diffy == y_grain)[1]]
  
  x_adj = 0 - (x_first - 0.5*x_grain - round((x_first - 0.5*x_grain)/x_grain)*x_grain) - 0.5*x_grain
  y_adj = 0 - (y_first - 0.5*y_grain - round((y_first - 0.5*y_grain)/y_grain)*y_grain) - 0.5*y_grain
  
  x_fuzzy_high = quadx + x_adj - round((quadx + x_adj)/x_grain)*x_grain + tol*x_grain
  x_fuzzy_low  = quadx + x_adj - round((quadx + x_adj)/x_grain)*x_grain - tol*x_grain
  x_in         = x_fuzzy_high >= 0 & x_fuzzy_low <= 0
  y_fuzzy_high = quady + y_adj - round((quady + y_adj)/y_grain)*y_grain + tol*y_grain
  y_fuzzy_low  = quady + y_adj - round((quady + y_adj)/y_grain)*y_grain - tol*y_grain
  y_in         = y_fuzzy_high >= 0 & y_fuzzy_low <= 0
  
  in_grid = x_in*y_in
  
  gridx = quadx[in_grid == 1]
  gridy = quady[in_grid == 1]
  
  newgridx = round((gridx + x_adj)/x_grain)*x_grain - x_adj
  newgridy = round((gridy + y_adj)/y_grain)*y_grain - y_adj
  
  outframe = data.frame(newgridx, newgridy, envframe[in_grid == 1, -c(xcol, ycol)])
  names(outframe)[1:2] = coord
  outframe
}


sampleQuad = function(env.grid, sp.scale, coord = c("X", "Y"), quadfilename = "Quad", tol = 0.01, writefile = TRUE)
{
	convert = FALSE
	env.grid = griddify(env.grid, tol = tol, coord = coord)
	if (any(lapply(env.grid, class) == "factor"))
	{
	  convert  = TRUE
	  out.grid = catConvert(env.grid)
	  env.grid = out.grid$X
	}
	x.col   = which(names(env.grid) == coord[1])
	y.col   = which(names(env.grid) == coord[2])
	res     = spatRes(env.grid, coord = coord)
	x.step  = res$x.step
	y.step  = res$y.step
	f.name = list()
	quad.list = list()
	for(i in 1:length(sp.scale))
	{
		i.scale = sp.scale[i]
		x.o = min(env.grid[,x.col]) - floor(min(env.grid[,x.col])/x.step)*x.step
		y.o = min(env.grid[,y.col]) - floor(min(env.grid[,y.col])/y.step)*y.step
		if (x.o/x.step > 0.5)
		{
			x.o = x.o - x.step
		}	
		if (y.o/y.step > 0.5)
		{
			y.o = y.o - y.step
		}

		is.on.scale   = abs((env.grid[,x.col]/i.scale) - round(env.grid[,x.col]/i.scale) - x.o/i.scale) + abs((env.grid[,y.col]/i.scale) - round(env.grid[,y.col]/i.scale) - y.o/i.scale) < 1.e-8
	  	dat.quad      = env.grid[is.on.scale,]
		### Get rid of machine error in coordinates
		dec.x = max(unlist(lapply(dat.quad[,x.col], decimalCount)))
		dec.y = max(unlist(lapply(dat.quad[,y.col], decimalCount)))
		dat.quad[,x.col] = unlist(lapply(dat.quad[, x.col], zapCoord, dec.x))
		dat.quad[,y.col] = unlist(lapply(dat.quad[, y.col], zapCoord, dec.y))
		#n.x           = floor((dat.quad[,x.col] - min(dat.quad[,x.col]))/x.step + 0.1)
		#n.y           = floor((dat.quad[,y.col] - min(dat.quad[,y.col]))/y.step + 0.1)
		#dat.quad[,x.col] = min(dat.quad[,x.col]) + n.x*x.step
		#dat.quad[,y.col] = min(dat.quad[,y.col]) + n.y*y.step

		if (is.na(quadfilename) == FALSE)
		{
			f.name[[i]] = paste(quadfilename, sp.scale[i], ".RData", sep = "")
			if (writefile == TRUE)
			{
			  save(dat.quad, file = f.name[[i]])
			  print(paste("Output saved in the file", f.name[[i]]))
			}
		}
		if (length(sp.scale) > 1)
		{
		  quad.list[[i]] = dat.quad
		}
	}
	if (convert == TRUE)
	{
		dat.quad = list(X = dat.quad, cat.names = out.grid$cat.names)
	}
	if (length(sp.scale) == 1)
		dat.quad
	else
		quad.list
}

ppmdat = function(sp.xy, sp.scale, back.xy, coord = c("X","Y"), sp.dat = getEnvVar(sp.xy = sp.xy, env.scale = sp.scale, env.grid = back.xy, coord = coord, writefile = writefile), sp.file = NA, quad.file = NA, datfilename = "PPMDat", writefile = TRUE)
{
	convert = FALSE
	if (is(sp.dat, "list"))
	{
		convert   = TRUE
		cat.names = sp.dat$cat.names
		sp.dat    = sp.dat$X
	}
	if (is.character(sp.xy) == TRUE)
	{
		sp.file = paste(sp.xy, "Env.RData")
		load(sp.file)
	}
	if (is.na(sp.file) == FALSE)
	{
		load(sp.file)
	}
	if (is.na(quad.file) == TRUE)
	{
		dat.quad = sampleQuad(env.grid = back.xy, sp.scale = sp.scale, coord = coord, writefile = writefile)
	}
	if (is.na(quad.file) != TRUE)
	{
		load(paste(quad.file, sp.scale, ".RData", sep = ""))
	}

	dat.quad$Pres = 0
	sp.dat$Pres   = 1
	
	quad.x.col = which(names(dat.quad) == coord[1])
	quad.y.col = which(names(dat.quad) == coord[2])
	sp.x.col   = which(names(sp.dat) == coord[1])
	sp.y.col   = which(names(sp.dat) == coord[2])
	sp.var     = setdiff(1:dim(sp.dat)[2], c(sp.x.col, sp.y.col))

	sp.data         = data.frame(sp.dat[,sp.x.col], sp.dat[,sp.y.col], sp.dat[,sp.var])
	names(sp.data)  = c("X", "Y", names(sp.dat)[sp.var])
	
	if (any(lapply(dat.quad, class) == "factor"))
	{
		dat.quad = catConvert(dat.quad)$X
	}
	add.var = setdiff(names(sp.dat), names(dat.quad))
	if (length(add.var) > 0)
	{
		quad.names = names(dat.quad)
		for (add in 1:length(add.var))
		{
			dat.quad = data.frame(dat.quad, 0)
		}
		names(dat.quad) = c(quad.names, add.var)

		dat.quad = dat.quad[,match(names(dat.quad), names(sp.dat))]
	}

	quad.var        = setdiff(1:dim(dat.quad)[2], c(quad.x.col, quad.y.col))
	quad.dat        = data.frame(dat.quad[,quad.x.col], dat.quad[,quad.y.col], dat.quad[,quad.var])
	names(quad.dat) = c("X", "Y", names(dat.quad)[quad.var])
	quad.dat        = quad.dat[,match(names(sp.dat), names(quad.dat))]
	dat.ppm         = rbind(sp.data, quad.dat)
	
	#dat.ppm$wt = ppmlasso:::weights(sp.data, quad.dat, coord)
	dat.ppm$wt = weights(sp.data, quad.dat, coord)

	dimnames(dat.ppm)[[1]] = 1:dim(dat.ppm)[1]
	if (convert == TRUE)
	{
		dat.ppm = list(X = dat.ppm, cat.names = cat.names)
	}
	if (is.na(datfilename) == FALSE)
	{
	  save.name = paste(datfilename, ".RData", sep = "")
	  if (writefile == TRUE)
	  {
	    save(dat.ppm, file = save.name)
	    print(paste("Output saved in the file", save.name))
	  }
	}
	dat.ppm
}

findRes = function(scales, lambda = 0, coord = c("X", "Y"), sp.xy, env.grid, formula, tol = 0.01, ...)
{
	form.1 = formula
	likelihoods = rep(NA, length(scales))
	#sp.data = env.var(sp.xy = sp.xy, env.scale = min(scales), env.grid = env.grid, coord = coord, file.name = "SpEnvData")
	sp.data = getEnvVar(sp.xy = sp.xy, env.scale = min(scales), env.grid = env.grid, coord = coord, envfilename = NA, tol = tol)
	for (sc in 1:length(scales))
	{
		formula = form.1
		if (lambda == 0)
		{
			data            = ppmdat(sp.xy = sp.xy, sp.scale = scales[sc], back.xy = env.grid, sp.dat = sp.data, sp.file = NA, quad.file = NA, writefile = FALSE)
			if (is(data, "list"))
			{
				use.form = as.character(formula)[2]
				cat.names = setdiff(unique(data$cat.names), NA)
				for (i in 1:length(cat.names))
				{
					use.form = gsub(cat.names[i], paste(names(data$X)[which(data$cat.names == cat.names[i])], collapse = " + "), use.form)
				}
				formula = as.formula(paste("~", use.form))
				data = data$X
			}
			glm.form = as.formula(paste("Pres/wt ~ ", as.character(formula)[2], sep = ""))
			glm.fit  = suppressWarnings(glm(glm.form, data = data, weights = data$wt, family = poisson()))
			eps   = 1.e-9
			while (glm.fit$deviance > glm.fit$null.dev)
			{
				glm.fit = suppressWarnings(glm(glm.form, data = data, weights = data$wt, family = poisson(), control = list(epsilon = eps)))
				eps     = eps/10
			}
			#glm.fit         = glm(glm.form, data = data, weights = data$wt, family = poisson(), control = list(epsilon = 1.e-12))
			likelihoods[sc] = sum(data$wt*(data$Pres/data$wt*log(glm.fit$fitted) - glm.fit$fitted)) - sum(log(1:sum(data$Pres > 0)))
		}
		if (lambda != 0)
		{
			sc.fit = ppmlasso(glm.form, sp.scale = scales[sc], lamb = lambda, data = ppmdat(sp.xy = sp.xy, sp.scale = scales[sc], back.xy = env.grid, sp.dat = sp.data, sp.file = NA, quad.file = NA, writefile = FALSE), ...)
			likelihoods[sc] = sc.fit$pen.likelihood[1]
		}
	}
	plot(scales, likelihoods, log = "x", type = "o", pch = 16, xlab = "Spatial Resolution", ylab = "Likelihood")
}

polyNames = function(X)
{
	coef.names = dimnames(X)[[2]]
	which.poly = grep("poly\\(", coef.names)
	if (length(which.poly) > 0)
	{
	  split1     = strsplit(coef.names[which.poly], "\\)")
	  #polyframe  =  plyr::ldply(split1, rbind)
	  #varframe   = plyr::ldply(strsplit(unlist(strsplit(as.character(polyframe[,1]), "poly\\("))[2*(1:(length(which.poly)))], ","), rbind)
	  polyframe  =  ldply(split1, rbind)
	  varframe   = ldply(strsplit(unlist(strsplit(as.character(polyframe[,1]), "poly\\("))[2*(1:(length(which.poly)))], ","), rbind)
	  varframe   = varframe[,-dim(varframe)[2]]
	  varframe   = data.frame(lapply(varframe, as.character), stringsAsFactors = FALSE)
	  if (dim(varframe)[1] == 1)
	  {
	    varframe = t(varframe)
	  }
	  #expframe   = plyr::ldply(strsplit(as.character(polyframe[,2]), "\\."), rbind)
	  expframe   = ldply(strsplit(as.character(polyframe[,2]), "\\."), rbind)
	  expframe   = data.frame(lapply(expframe, as.character), stringsAsFactors = FALSE)
	  expframe[is.na(expframe)] = 0
	  vframe = varframe[,1:dim(expframe)[2]]
	  nameframe  = matrix(paste(as.matrix(vframe), "_EXP", as.matrix(expframe), sep = ""), nrow(vframe), ncol(vframe))
	  nameframe[expframe == "0"] = ""
	  nameframe[expframe == "1"] = as.character(vframe[expframe == 1])
	  nameframe = gsub(" ", "", nameframe)
	  nameframe  = as.data.frame(nameframe)
	  if (is.null(dim(nameframe)) == FALSE)
	  {
	    names.out = apply(nameframe, 1, function(row) paste(row[nzchar(row)], collapse = "_x_"))
	  }
	  if (is.null(dim(nameframe)) == TRUE)
	  {
	    names.out = nameframe
	  }
	  id.out    = which.poly
	}
	else
	{
	  names.out = NULL
	  id.out = NULL
	}
	return(list(names = names.out, ids = id.out))
}

ppmlasso = function(formula, sp.xy, env.grid, sp.scale, coord = c("X", "Y"), data = ppmdat(sp.xy = sp.xy, sp.scale = sp.scale, back.xy = env.grid, coord = c("X","Y"), sp.file = NA, quad.file = NA, datfilename = "PPMDat", writefile = writefile), lamb = NA, n.fits = 200, ob.wt = NA, criterion = "bic", alpha = 1, family = "poisson", tol = 1.e-9, gamma = 0, init.coef = NA, mu.min = 1.e-16, mu.max = 1/mu.min, r = NA, interactions = NA, availability = NA, max.it = 25, min.lamb = -10, standardise = TRUE, n.blocks = NA, block.size = sp.scale*100, seed = 1, writefile = TRUE)
{
	error.flag = FALSE
	formula.out = formula
	if (is(data, "list"))
	{
		use.form = as.character(formula)[2]
		cat.names = setdiff(unique(data$cat.names), NA)
		for (i in 1:length(cat.names))
		{
			use.form = gsub(cat.names[i], paste(names(data$X)[which(data$cat.names == cat.names[i])], collapse = " + "), use.form)
		}
		formula = as.formula(paste("~", use.form))
		data = data$X
	}
	wt.calc = FALSE
	if (is.na(ob.wt) == TRUE)
	{
		wt.calc = TRUE
		ob.wt   = data$wt
	}

	call = match.call()
	mf   = model.frame(formula, data = data)
	mt   = attr(mf, "terms")
   	y    = data$Pres/data$wt
	X    = if (!is.empty.model(mt)) 
        model.matrix(mt, mf)
   	 else matrix(, NROW(y), 0L)
	
	if (any(X[,1] != 1))
	{
		X = as.matrix(cbind(1, X))
	}
	cut.var = which(apply(X, 2, max) == 0)
	if (length(cut.var) > 0)
	{
		X       = X[, -cut.var]
	}
	Xnames  = polyNames(X)
	#Xnames = colnames(X)
	if (is.null(Xnames) == FALSE)
	{
		dimnames(X)[[2]][Xnames$ids] = Xnames$names
	}
	dimnames(X)[[2]][1] = "Intercept"

	area.int = FALSE
	raw.int  = NA
	if (family == "area.inter")
	{
		family   = "poisson"
		area.int = TRUE
		if (is.na(interactions) == TRUE)
		{
			interactions = pointInteractions(data, r, availability)
		}
		raw.int = interactions
	}

	s.means = NULL
	s.sds   = NULL

	if (standardise == TRUE)
	{
		stand.X = standardiseX(X[,-1])
		X       = as.matrix(cbind(1, stand.X$X))
		dimnames(X)[[2]][1] = "Intercept"
		s.means = stand.X$dat.means
		s.sds   = stand.X$dat.sds
		if (area.int == TRUE)
		{
			interactions = standardiseX(interactions)$X
		}
	}

	if (family == "poisson")
	{
        	family = get(family, mode = "function", envir = parent.frame())
	}
   	if (is.function(family))
	{
        	family = family()
	}
	if (family$family == "binomial")
	{
		mu.max  = 1 - mu.min
	}
	score.0    = scoreIntercept(y, X, ob.wt = ob.wt, area.int = area.int, int = interactions, family)

	if (is.na(init.coef))
	{
		gamma = 0
	}
	if (gamma == 0)
	{
		adapt.weights = rep(1, dim(X)[2])
	}
	if (gamma != 0)
	{
		adapt.weights = 1/abs(init.coef)^gamma
	}

	cv = rep(0, dim(data)[1])

	if (criterion == "blockCV")
	{
		cv      = blocks(n.blocks, block.size, data, seed = seed)
		pred.mu = matrix(NA, dim(data)[1], n.fits)
	}
	data.all = data
	y.all    = y
	X.all    = X
	wt.all   = ob.wt
	if (area.int == TRUE)
	{
		interactions.all = interactions
	}

	for (cv.i in 1:length(unique(cv)))
	{
	dat.test = data[cv == cv.i,]
	data     = data[cv != cv.i,]
	data$wt  = weights(data[data$Pres == 1,], data[data$Pres == 0,], coord)
	y        = data$Pres/data$wt
	X        = X[cv != cv.i,]
	ob.wt    = ob.wt[cv != cv.i]
	if (wt.calc == TRUE)
	{
		ob.wt = data$wt
	}
	if (area.int == TRUE)
	{
		interactions = interactions[cv != cv.i]
	}

	if (is.na(lamb) == TRUE)
	{
		new.score  = abs(score.0/adapt.weights)
		sub.score  = new.score[-1]
		max.lambda = max(sub.score[is.infinite(sub.score) == FALSE])
		if (is.na(max.lambda) == FALSE)
		{
			lambs = sort(exp(seq(min.lamb, log(max.lambda + 1.e-5), length.out = n.fits)), decreasing = TRUE)
		}
	}
	if (is.na(lamb) == FALSE)
	{
		lambs = lamb
	}

	if (criterion == "blockCV")
	{
		lambda.max = max(abs(scoreIntercept(data$Pres/data$wt, X, ob.wt = data$wt, area.int = area.int, int = interactions, family = poisson())[-1]))	
		lambs      = exp(seq(0, -12, length.out = n.fits))*lambda.max
	}

	n.pres     = sum(y > 1.e-8)
	n.var      = dim(X)[2] - 1
	
	if (criterion == "msi")
	{
		lambs = max.lambda/sqrt(n.pres)
	}

	if (area.int == TRUE)
	{
		X.0 = as.matrix(cbind(X, interactions))
	}
	if (area.int != TRUE)
	{
		X.0 = X
	}

	mod.0     = glm(y ~ X.0[,-1], family = family, weights = ob.wt)
	#coefs     = matrix(NA, (dim(X)[2] + area.int), (length(lambs) + 1), dimnames = list(dimnames(X)[[2]]))
	coefs     = matrix(NA, (dim(X)[2]), (length(lambs) + 1), dimnames = list(dimnames(X)[[2]]))
	if (area.int == 1)
	{
		coefs     = matrix(NA, (dim(X)[2] + area.int), (length(lambs) + 1), dimnames = list(c(dimnames(X)[[2]], "Interaction")))
	}
	num.param = rep(NA, (length(lambs) + 1))
	gcvs      = rep(NA, (length(lambs) + 1))
	aics      = rep(NA, (length(lambs) + 1))
	hqcs      = rep(NA, (length(lambs) + 1))
	bics      = rep(NA, (length(lambs) + 1))
	devs      = rep(NA, (length(lambs) + 1))
	ll        = rep(NA, (length(lambs) + 1))
	pll       = rep(NA, (length(lambs) + 1))
	nlgcvs    = rep(NA, (length(lambs) + 1))
	offset    = c(log(mean(y)), rep(0, n.var), rep(1, area.int))

	if (any(is.na(adapt.weights) == FALSE))
	{
		if (sum(is.infinite(adapt.weights[-1])) != (length(adapt.weights) - 1 - area.int))
		{
			it.max = 100
			for (reg.path in 1:length(lambs))
			{
				mod  = try(singleLasso(y, X, lamb = lambs[reg.path], ob.wt = ob.wt, alpha = alpha, b.init = offset, family = family, tol = 1.e-9, gamma = gamma, init.coef = init.coef, area.int = area.int, interactions = interactions, max.it = it.max, standardise = FALSE), TRUE)
				if(is(mod, "try-error"))
				{
					break
				}
				if (any(is.na(mod$b)))
				{
					break
				}
				coefs[,reg.path]    = mod$b
				gcvs[reg.path]      = mod$GCV
				aics[reg.path]      = mod$AIC
				hqcs[reg.path]      = mod$HQC
				bics[reg.path]      = mod$BIC
				devs[reg.path]      = mod$dev
				ll[reg.path]        = mod$like
				pll[reg.path]       = mod$pen
				num.param[reg.path] = sum(abs(mod$b) > 1.e-7) - 1 - area.int
				offset = mod$b
				it.max = max.it
				cat(paste("Fitting Models:", reg.path, "of", length(lambs), "\n"))
				flush.console()
			}

			num.done = length(lambs) + 1 - sum(is.na(aics))
			s.denom  = sum(abs(mod.0$coef[-c(1, (n.var + 2))]))
	
			if (n.var > n.pres)
			{
				s.denom = sum(abs(coefs[-c(1, (n.var + 2)), num.done]))
			}

			if (any(is.na(mod.0$coef)) == TRUE)
			{
				s.denom = sum(abs(coefs[-c(1, (n.var + 2)), num.done]))
			}

			ss     = rep(NA, (length(lambs) + 1))
			nlgcvs = rep(NA, (length(lambs) + 1))

			for (i.nlgcv in 1:num.done)
			{
				s.num       = sum(abs(coefs[-c(1, (n.var + 2)), i.nlgcv]))
				s           = s.num/s.denom
				ss[i.nlgcv] = s
				nlgcv       = devs[i.nlgcv]/(n.pres*(1 - (area.int + n.var*s)/n.pres)^2)
				if (area.int + n.var*s > n.pres)
				{
					nlgcv = NA
				}
				nlgcvs[i.nlgcv] = nlgcv
			}
		}
		
		if (sum(is.infinite(adapt.weights[-1])) == (length(adapt.weights) - 1 - area.int))
		{
			coefs     = matrix(rep(init.coef, (length(lambs) + 1)), length(adapt.weights), (length(lambs) + 1))
			gcvs      = rep(1, (length(lambs) + 1))
			aics      = rep(1, (length(lambs) + 1))
			bics      = rep(1, (length(lambs) + 1))
			hqcs      = rep(1, (length(lambs) + 1))
			nlgcvs    = rep(1, (length(lambs) + 1))
			devs      = rep(1, (length(lambs) + 1))
			ll        = rep(1, (length(lambs) + 1))
			num.param = rep(0, (length(lambs) + 1))
			num.done  = length(lambs) + 1 - sum(is.na(aics))
		}

		criterion.matrix        = data.frame(aics, bics, hqcs, gcvs, nlgcvs)
		names(criterion.matrix) = c("AIC", "BIC", "HQC", "GCV", "NLGCV")

		lambs[(length(lambs) + 1)] = 0
		coefs[,length(lambs)] = mod.0$coef
		if (gamma != 0)
		{
			coefs[,length(lambs)] = init.coef
		}
		num.param[(length(lambs) + 1)] = length(adapt.weights) - sum(is.infinite(adapt.weights[-c(1, (n.var + 2))])) - 1 - area.int
		
		meth.id  = paste(criterion, "s", sep = "")

		if (criterion == "msi" | criterion == "blockCV")
		{
			choice.id = 1
		}
		if (criterion != "msi" & criterion != "blockCV")
		{
			choice.id = max(which.min(get(meth.id)))
		}

		lambda.hat = lambs[choice.id]
		beta.hat   = coefs[,choice.id]
		eta.hat    = X.0 %*% beta.hat
		mu.hat     = family$linkinv(eta.hat)
		like.hat   = ll[choice.id]

		assign(paste("coefs.", cv.i, sep = ""), coefs)
		

		if (criterion == "blockCV")
		{
			for (i in 1:length(lambs) - 1)
			{ #Fix this for predicting to test locations
				fam.fit = poisson()
				if (area.int == TRUE)
				{
					fam.fit = "area.inter"
				}
				cv.fit = list(beta = coefs[,i], s.means = s.means, s.sds = s.sds, family = fam.fit, pt.interactions = interactions, formula = formula, mu = rep(0, dim(data)[1]))
				if (area.int == TRUE)
				{
					pred.mu[cv == cv.i, i] = predict.ppmlasso(cv.fit, newdata = dat.test, interactions = interactions.all[cv == cv.i])
				}
				if (area.int != TRUE)
				{
					pred.mu[cv == cv.i, i] = predict.ppmlasso(cv.fit, newdata = dat.test)
				}
			}
		}
	}

	data  = data.all
	y     = y.all
	X     = X.all
	ob.wt = wt.all
	if (area.int == TRUE)
	{
		interactions = interactions.all
	}

	} # cv.i

	if (criterion == "blockCV")
	{
		l.vec       = apply(pred.mu, 2, calculateUnpenalisedLikelihood, ob.wt = data$wt, y = data$Pres/data$wt)
		coef.mat    = matrix(NA, n.blocks, dim(coefs)[1])
		for (i in 1:n.blocks)
		{
			coefs.i       = get(paste("coefs.", i, sep = ""))
			coef.mat[i, ] = coefs.i[, which.max(l.vec)]
		}
		coef.av = apply(coef.mat, 2, mean, na.rm = TRUE)
		lambda.mult = exp(seq(0, -12, length.out = n.fits))[which.max(l.vec)]
		if (area.int != TRUE)
		{
			lambda.max = max(abs(scoreIntercept(data$Pres/data$wt, X, ob.wt = data$wt, family = poisson())[-1]))
		}	

		if (area.int == TRUE)
		{
			lambda.max = max(abs(scoreIntercept(data$Pres/data$wt, X, ob.wt = data$wt, family = poisson(), area.int = TRUE, int = scale(interactions))[-1]))
		}
		final.fit = try(singleLasso(y.all, X.all, lamb = lambda.mult*lambda.max, ob.wt = wt.all, alpha = alpha, b.init = coef.av, family = family, tol = 1.e-9, gamma = gamma, init.coef = init.coef, area.int = area.int, interactions = interactions, max.it = it.max, standardise = FALSE), TRUE)
		coefs     = final.fit$b
		lambs     = lambda.mult*lambda.max
		gcvs      = NA
		aics      = NA
		hqcs      = NA
		bics      = NA
		nlgcvs    = NA
		devs      = final.fit$dev
		ll        = final.fit$like
		pll       = final.fit$pen
		num.param = sum(abs(final.fit$b) > 1.e-7) - 1 - area.int

		if (area.int == TRUE)
		{
			X.0 = as.matrix(cbind(X, interactions))
		}
		if (area.int != TRUE)
		{
			X.0 = X
		}
		
		lambda.hat = lambs
		beta.hat   = coefs
		eta.hat    = X.0 %*% beta.hat
		mu.hat     = family$linkinv(eta.hat)
		like.hat   = ll
		criterion.matrix = data.frame(aics, bics, hqcs, gcvs, nlgcvs)
	}

	family.out = family$family
	if (area.int == TRUE)
	{	
		family.out = "area.inter"
	}
	output = list(betas = coefs, lambdas = lambs, likelihoods = ll, pen.likelihoods = pll, lambda = lambda.hat, beta = beta.hat, mu = mu.hat, likelihood = like.hat, criterion = criterion, family = family.out, gamma = gamma, alpha = alpha, init.coef = init.coef, criterion.matrix = criterion.matrix, data = X.0, pt.interactions = raw.int, wt = ob.wt, pres = data$Pres, x = data$X, y = data$Y, r = r, call = call, formula = formula.out, s.means = s.means, s.sds = s.sds, cv.group = cv, n.blocks = n.blocks)
	class(output) = c("ppmlasso", "list")
	output
}

interp = function(sp.xy, sp.scale, f, back.xy, coord = c("X","Y"))
{
	options(scipen = 999)
	x.dat   = sp.xy[,which(names(sp.xy) == coord[1])]
	y.dat   = sp.xy[,which(names(sp.xy) == coord[2])]
	x.back  = back.xy[,which(names(back.xy) == coord[1])]
	y.back  = back.xy[,which(names(back.xy) == coord[2])]
	res     = spatRes(back.xy, coord = coord)

	grid    = data.table(x.back, y.back, f, key = c("x.back", "y.back"))

	ux    = sort(unique(x.back))
	uy    = sort(unique(y.back))

	x.col   = which(names(back.xy) == coord[1])
	y.col   = which(names(back.xy) == coord[2])

	x.step  = res$x.step
	y.step  = res$y.step

	x.o = min(back.xy[,x.col]) - floor(min(back.xy[,x.col])/x.step)*x.step
	y.o = min(back.xy[,y.col]) - floor(min(back.xy[,y.col])/y.step)*y.step

	x.1   = floor((x.dat - x.o)/sp.scale)*sp.scale + x.o
	y.1   = floor((y.dat - y.o)/sp.scale)*sp.scale + y.o
	x.2   = pmin(x.1 + sp.scale, max(ux))
	y.2   = pmin(y.1 + sp.scale, max(uy))

	w11   = (x.2 - x.dat)*(y.2 - y.dat)/((x.2 - x.1)*(y.2 - y.1))
	w12   = (x.2 - x.dat)*(y.dat - y.1)/((x.2 - x.1)*(y.2 - y.1))
	w21   = (x.dat - x.1)*(y.2 - y.dat)/((x.2 - x.1)*(y.2 - y.1))
	w22   = (x.dat - x.1)*(y.dat - y.1)/((x.2 - x.1)*(y.2 - y.1))

	f11   = grid[list(x.1, y.1)]$f
	f12   = grid[list(x.1, y.2)]$f
	f21   = grid[list(x.2, y.1)]$f
	f22   = grid[list(x.2, y.2)]$f

	c11 = 1 - is.na(f11)
	c12 = 1 - is.na(f12)
	c21 = 1 - is.na(f21)
	c22 = 1 - is.na(f22)

	env.wt.mat = cbind(f11*w11*c11, f12*w12*c12, f21*w21*c21, f22*w22*c22) 
	
	f.interp = apply(env.wt.mat, 1, sum, na.rm = TRUE)/(w11*c11 + w12*c12 + w21*c21 + w22*c22)
	f.interp
}

getEnvVar = function(sp.xy, env.grid, env.scale, coord = c("X", "Y"), envfilename = "SpEnvData", tol = 0.01, writefile = TRUE)
{
  env.grid = griddify(env.grid, tol = tol, coord = coord)
	convert = FALSE
	if (any(lapply(env.grid, class) == "factor"))
	{
		convert  = TRUE
		out.grid = catConvert(env.grid)
		env.grid = out.grid$X
	}
	x.dat   = sp.xy[,which(names(sp.xy) == coord[1])]
	y.dat   = sp.xy[,which(names(sp.xy) == coord[2])]
	x.back  = env.grid[,which(names(env.grid) == coord[1])]
	y.back  = env.grid[,which(names(env.grid) == coord[2])]
	x.col   = which(names(env.grid) == coord[1])
	y.col   = which(names(env.grid) == coord[2])
	var.col = setdiff(1:dim(env.grid)[2], c(x.col, y.col))
	s.res   = spatRes(env.grid)
	
	sp.dat        = as.data.frame(matrix(NA, length(x.dat), length(var.col)))
	names(sp.dat) = names(env.grid[var.col])

	for (var in 1:length(var.col))
	{
		loop.scale = min(c(s.res$x.step, s.res$y.step))
		loc        = which(is.na(sp.dat[,var]))
		while(sum(is.na(sp.dat[,var])) > 0)
		{
			loc = which(is.na(sp.dat[,var]))
			sp.dat[loc, var] = interp(sp.xy[loc,], loop.scale, env.grid[,var.col[var]], env.grid, coord = c("X","Y"))
			loop.scale = loop.scale*2
		}
		cat(paste("Calculating species environmental data for variable:", names(sp.dat)[var], "\n"))
		flush.console()
	}
	
	sp.dat = data.frame(x.dat, y.dat, sp.dat)
	names(sp.dat)[1:2] = c("X", "Y")
	if (is.na(envfilename) == FALSE)
	{
		save.name = paste(envfilename, ".RData", sep = "")
		if (writefile == TRUE)
		{
		  save(sp.dat, file = save.name)
		  print(paste("Output saved in the file", save.name))
		}
	}
	if (convert == TRUE)
	{
		sp.dat = list(X = sp.dat, cat.names = out.grid$cat.names)
	}
	sp.dat
}

makeMask = function(dat.ppm)
{
	if (is(dat.ppm, "list"))
	{
		dat.ppm = dat.ppm$X
	}
	dat.q = dat.ppm[dat.ppm$Pres == 0,]
	ux = sort(unique(dat.q$X))
	uy = sort(unique(dat.q$Y))
	nx = length(ux)
	ny = length(uy)

	col.ref = match(dat.q$X, ux)
	row.ref = match(dat.q$Y, uy)

	all.vec          = rep(0, max(row.ref)*max(col.ref))
	vec.ref          = (col.ref - 1)*max(row.ref) + row.ref
	all.vec[vec.ref] = 1
	mask.out         = matrix(all.vec, max(row.ref), max(col.ref), dimnames = list(uy, ux))
	mask.out
}


pointInteractions = function(dat.ppm, r, availability = NA)
{
	if (is(dat.ppm, "list"))
	{
		dat.ppm = dat.ppm$X
	}
	if (any(is.na(availability)))
	{
		#grain = min(0.5, r/2)
		#min.x = min(dat.ppm$X) - r - grain
		#max.x = max(dat.ppm$X) + r + grain
		#min.y = min(dat.ppm$Y) - r - grain
		#max.y = max(dat.ppm$Y) + r + grain

		#availability = matrix(1, 1 + (max.y - min.y)/grain, 1 + (max.x - min.x)/grain)
		#rownames(availability) = seq(min.y, max.y, grain)
		#colnames(availability) = seq(min.x, max.x, grain)
		availability = makeMask(dat.ppm)
	}	

	cat(paste("Calculating point interactions", "\n"))
	flush.console()
	occupied = matrix(0, dim(availability)[1], dim(availability)[2])
	rownames(occupied) = rownames(availability)
	colnames(occupied) = colnames(availability)

	x.mat = availability
	y.mat = availability

	for (i in 1:dim(x.mat)[1])
	{
		x.mat[i,] = as.numeric(colnames(availability))
	}
	for (i in 1:dim(y.mat)[2])
	{
		y.mat[,i] = as.numeric(rownames(availability))
	}

	grain = as.numeric(colnames(availability)[2]) - as.numeric(colnames(availability)[1])

	pres.x = dat.ppm$X[dat.ppm$Pres > 0]
	pres.y = dat.ppm$Y[dat.ppm$Pres > 0]

	quad.x = dat.ppm$X[dat.ppm$Pres == 0]
	quad.y = dat.ppm$Y[dat.ppm$Pres == 0]

	quad.int = rep(0, length(quad.x))

	for (i in 1:length(pres.x))
	{
		sub.col = which(as.numeric(colnames(occupied)) >= pres.x[i] - (r + grain) & as.numeric(colnames(occupied)) <= pres.x[i] + (r + grain))
		sub.row = which(as.numeric(rownames(occupied)) >= pres.y[i] - (r + grain) & as.numeric(rownames(occupied)) <= pres.y[i] + (r + grain))

		sub.occ = occupied[sub.row, sub.col]
		sub.x   = x.mat[sub.row, sub.col]
		sub.y   = y.mat[sub.row, sub.col]

		sub.occ[(sub.x - pres.x[i])^2 + (sub.y - pres.y[i])^2 < r^2] = sub.occ[(sub.x - pres.x[i])^2 + (sub.y - pres.y[i])^2 < r^2] + 1
		occupied[sub.row, sub.col] = sub.occ
		quad.cells           = which((quad.x - pres.x[i])^2 + (quad.y - pres.y[i])^2 <= (2*r)^2)
		quad.int[quad.cells] = quad.int[quad.cells] + 1
	}

	int.q = rep(0, length(quad.int))

	for (quad.i in which(quad.int > 0))
	{
		sub.col = which(as.numeric(colnames(occupied)) >= quad.x[quad.i] - (r + grain) & as.numeric(colnames(occupied)) <= quad.x[quad.i] + (r + grain))
		sub.row = which(as.numeric(rownames(occupied)) >= quad.y[quad.i] - (r + grain) & as.numeric(rownames(occupied)) <= quad.y[quad.i] + (r + grain))

		sub.occ  = occupied[sub.row, sub.col]
		sub.availability = availability[sub.row, sub.col]
		sub.x    = x.mat[sub.row, sub.col]
		sub.y    = y.mat[sub.row, sub.col]

		sub.cell = (sub.x - quad.x[quad.i])^2 + (sub.y - quad.y[quad.i])^2 <= r^2 & sub.availability > 0

		#int.q[quad.i] = sum(sub.occ[sub.cell] > 0)/sum(sub.cell)
		int.q[quad.i] = sum(sub.occ[sub.cell] > 0, na.rm = TRUE)/sum(sub.cell, na.rm = TRUE)
	}

	int.p = rep(0, length(pres.x))

	for (pres.i in 1:length(pres.x))
	{
		sub.col = which(as.numeric(colnames(occupied)) >= pres.x[pres.i] - (r + grain) & as.numeric(colnames(occupied)) <= pres.x[pres.i] + (r + grain))
		sub.row = which(as.numeric(rownames(occupied)) >= pres.y[pres.i] - (r + grain) & as.numeric(rownames(occupied)) <= pres.y[pres.i] + (r + grain))

		sub.occ  = occupied[sub.row, sub.col]
		sub.availability = availability[sub.row, sub.col]
		sub.x    = x.mat[sub.row, sub.col]
		sub.y    = y.mat[sub.row, sub.col]

		sub.cell = (sub.x - pres.x[pres.i])^2 + (sub.y - pres.y[pres.i])^2 <= r^2 & sub.availability > 0

		#int.p[pres.i] = sum(sub.occ[sub.cell] > 1)/sum(sub.cell)
		int.p[pres.i] = sum(sub.occ[sub.cell] > 1, na.rm = TRUE)/sum(sub.cell, na.rm = TRUE)
	}

	interactions = c(int.p, int.q)
	interactions
}

plotPath = function(fit, colors = c("gold", "green3", "blue", "brown", "pink"), logX = TRUE)
{
	best.models = apply(fit$criterion.matrix, 2, which.min) #See which fit optimises other criteria
	unique.beta = unique(best.models)
	names = rep(NA, length(unique.beta))
	for (i in 1:length(unique.beta))
	{
		names[i] = paste(names(best.models[best.models == unique.beta[i]]), collapse = "/")
	}

	min.y = min(apply(fit$betas, 1, min))
	max.y = max(apply(fit$betas, 1, max))
	
	keep_models = if (logX == TRUE) 1:(length(fit$lambdas) - 1) else 1:length(fit$lambdas)
	
	if (logX == TRUE)
	{
	  plot(fit$lambdas[keep_models], fit$betas[1, keep_models], log = "x", type = "l", ylim = c(min.y, max.y), xlab = "LASSO penalty", ylab = "Coefficients") 
	}
	else
	{
	  plot(fit$lambdas[keep_models], fit$betas[1, keep_models], type = "l", ylim = c(min.y, max.y), xlab = "LASSO penalty", ylab = "Coefficients") 
	}
	  
	for (i in 2:dim(fit$betas)[1])
	{
		points(fit$lambdas[keep_models], fit$betas[i, keep_models], type = "l")
	}

	for (i in 1:length(unique.beta))
	{
		abline(v = fit$lambdas[unique.beta[i]], lwd = 3, col = colors[i])
	}	

	legend("topright", names, lwd = rep(3, length(unique.beta)), col = colors[1:length(unique.beta)])
}

plotFit = function(fit, pred.data = data.frame(X = fit$x[fit$pres == 0], Y = fit$y[fit$pres == 0], 1, scale(fit$data[fit$pres == 0, -1], center = -fit$s.means/fit$s.sds, scale = 1/fit$s.sds)), 
	coord = c("X", "Y"), asp = "iso", ylab = "", xlab = "", col.regions = heat.colors(1024)[900:1], 
	cuts = length(col.regions), cex = 1.4, main.text = paste(toupper(fit$criterion), "fit"), cex.color = 1.4)
{
	x.col    = which(names(pred.data) == coord[1])
	y.col    = which(names(pred.data) == coord[2])
	pred.int = predict(fit, newdata = pred.data[,-c(x.col, y.col)])
	levelplot(pred.int ~ pred.data[,x.col] + pred.data[,y.col], asp = asp, 
		ylab = ylab, xlab = xlab, col.regions = col.regions, cuts = cuts, 
		main = list(main.text, cex = cex), 
		scales = list(y = list(draw = FALSE), x = list(draw = FALSE)), 
		cex = cex, colorkey = list(labels = list(cex = cex.color)))

}

catConvert = function(env.frame)
{
  classes   = lapply(env.frame, class)
  is.cat    = which(classes == "factor")
  cat.names = names(env.frame)[is.cat]
  frame.out = env.frame[, which(classes != "factor")]
  cat.out   = rep(NA, dim(frame.out)[2])
  for (v in 1:length(is.cat))
  {
    current.dim = dim(frame.out)[2]
    factors = sort(unique(env.frame[, is.cat[v]]))
    fac.frame = c()
    for (i in 1:length(factors))
    {
      fac.frame = cbind(fac.frame, as.numeric(env.frame[, is.cat[v]] == factors[i]))
    }
    fac.frame = data.frame(fac.frame)
    frame.out = data.frame(frame.out, fac.frame)
    names(frame.out)[(current.dim + 1):dim(frame.out)[2]] = paste(cat.names[v], ".", factors, sep = "")
    cat.out = c(cat.out, rep(cat.names[v], dim(fac.frame)[2]))
  }
  return(list(X = frame.out, cat.names = cat.out))
}

catFrame = function(cat.mat)
{
  cat.combine = setdiff(unique(cat.mat$cat.names), NA)
  frame.out   = cat.mat$X[, is.na(cat.mat$cat.names)]
  for (i in 1:length(cat.combine))
  {
    frame.cat = cat.mat$X[, which(cat.mat$cat.names == cat.combine[i])]
    frame.out = within(frame.out, {assign(cat.combine[i], frame.cat)})
  }
  frame.out
}

standardiseX = function(mat)
{
  X = scale(as.matrix(mat))
  dat.means = apply(as.matrix(mat), 2, mean, na.rm = TRUE)
  dat.sds   = apply(as.matrix(mat), 2, sd, na.rm = TRUE)
  return(list(X = X, dat.means = dat.means, dat.sds = dat.sds))
}

muFromEta = function(eta, family, mu.min = 1.e-16, mu.max = 1/mu.min)
{
  mu = family$linkinv(eta)
  mu[mu < mu.min] = mu.min
  mu[mu > mu.max] = mu.max
  mu
}	

etaFromMu = function(mu, family, mu.min = 1.e-16, mu.max = 1/mu.min, eta.min, eta.max)
{
  eta = family$linkfun(mu)
  eta[eta < eta.min] = eta.min
  eta[eta > eta.max] = eta.max
  eta
}

irlsUpdate = function(y, X, ob.wt, is.in, signs, eta, mu, alpha, lambda, beta.old, penalty = FALSE, family, mu.min = 1.e-16, mu.max = 1/mu.min, eta.min, eta.max, tol = tol)
{
  if (penalty == FALSE)
  {
    is.in  = rep(TRUE, dim(X)[2])
    lambda = rep(0, dim(X)[2])
    signs  = rep(0, dim(X)[2])
  }
  
  vari  = family$variance(mu)
  deriv = 1/vari
  z     = eta + (y - mu)*deriv
  weii  = ob.wt*1/(deriv^2*vari)
  
  Xw       = t(as.vector(weii) * t(t(X[,is.in])))
  Xw.s     = t(as.vector(weii) * t(t(X)))
  epsilon  = diag(rep(sqrt(tol), sum(is.in)))
  if (sum(is.in) == 1)
  {
    epsilon = sqrt(tol)
  }
  
  xwx = Xw %*% X[,is.in] + epsilon
  qx = qr(xwx)
  dim.X = dim(xwx)
  if (qx$rank >= dim.X[1])
  {
    if (dim.X[1] == 1)
    {
      beta.new = solve(xwx + (1 - alpha)*lambda[is.in]) %*% (Xw %*% z - as.matrix(alpha * lambda[is.in] * signs[is.in]))
    }
    if (dim.X[1] > 1)
    {
      beta.new = solve(xwx + diag((1 - alpha)*lambda[is.in])) %*% (Xw %*% z - as.matrix(alpha * lambda[is.in] * signs[is.in]))
    }
    
    eta.new  = as.matrix(X[,is.in]) %*% as.matrix(beta.new)
    eta.new[eta.new < eta.min] = eta.min
    eta.new[eta.new > eta.max] = eta.max
    
    mu.new = family$linkinv(eta.new)
    mu.new[mu.new < mu.min] = mu.min
    mu.new[mu.new > mu.max] = mu.max
    return(list(mu = mu.new, beta = beta.new, eta = eta.new, wt = Xw, XwX = xwx, s.wt = Xw.s, deriv = deriv, v = vari, error = "None"))
  }
  if (qx$rank < dim.X)
  {
    error.flag = "Singular matrix"
    return(list(error = error.flag))
  }
}

calculateLikelihood = function(X, family, ob.wt, mu, y, alpha, lambda, beta, penalty = FALSE)
{
  if (penalty == FALSE)
  {
    lambda = rep(0, dim(X)[2])
  }
  if (family$family == "poisson")
  {
    like = sum(ob.wt*(y*log(mu) - mu)) - sum(log(1:sum(y > 0))) - sum(alpha * as.vector(lambda)*abs(beta)) - sum(0.5 * (1 - alpha) * as.vector(lambda)*beta^2)
  }
  if (family$family == "binomial")
  {
    like = sum(ob.wt*(y*log(mu) + (1 - y)*log(1 - mu))) - sum(as.vector(lambda)*abs(beta))
  }
  if (family$family == "gaussian")
  {
    like = sum(ob.wt*(y*mu - 0.5*mu^2 - 0.5*y^2 - 0.5*log(2*pi))) - sum(as.vector(lambda)*abs(beta))
  }
  like
}

calculateGCV = function(y, X, ob.wt, b.lasso, lambda, alpha = alpha, unp.likelihood = unp.likelihood, penalty = TRUE, family = "poisson", mu.min = 1.e-16, mu.max = 1.e16, eta.min = log(1.e-16), eta.max = log(1.e16), tol = 1.e-9, area.int = FALSE)
{
  is.in        = abs(b.lasso) > 1.e-7
  signs        = sign(b.lasso)
  eta          = X %*% b.lasso
  mu           = exp(eta)
  n.pres       = sum(y > 0)
  wt           = irlsUpdate(y, X, ob.wt, is.in, signs, eta, mu, alpha, lambda, b.lasso, penalty = TRUE, family = family, mu.min, mu.max, eta.min, eta.max, tol)$wt
  xwx          = wt %*% X[,is.in]
  bpinv        = rep(0, length(b.lasso))
  bpinv[is.in] = 1/abs(b.lasso[is.in])
  ginv         = diag(bpinv)	
  
  if (sum(is.in) > 1)
  {
    eff.df = sum(diag(solve(xwx + diag(as.vector(lambda[is.in])) %*% ginv[is.in, is.in]) %*% xwx)) + area.int
  }
  if (sum(is.in) == 1)
  {
    eff.df = sum(diag(solve(xwx + diag(as.matrix(lambda[is.in])) %*% ginv[is.in, is.in]) %*% xwx)) + area.int
  }
  
  if (family$family == "poisson")
  {
    dev = -2*(sum(ob.wt*(y*log(mu) - mu)) + sum(y > 0)) 
  }
  
  if (family$family == "binomial")
  {
    dev = -2*unp.likelihood
  }
  
  if (family$family == "gaussian")
  {
    dev = -2*(unp.likelihood - calculateLikelihood(X, family, ob.wt[y > 0], y[y > 0], y, alpha, lambda, beta = b.lasso, penalty = FALSE))
  }
  
  if (eff.df < n.pres)
  {
    gcv = dev/(n.pres*(1 - (eff.df)/n.pres)^2)
  }
  if (eff.df >= n.pres)
  {
    gcv = NA
  }
  return(list(gcv = gcv, dev = dev, eff.df = eff.df))
}

scoreIntercept = function(y, X.des, ob.wt = rep(1, length(y)), area.int = FALSE, int = NA, family)
{
  if (area.int == FALSE)
  {
    int.mod = glm(y ~ 1, family = family, weights = ob.wt)
  }
  if (area.int != FALSE)
  {
    int.mod = glm(y ~ 1 + int, family = family, weights = ob.wt)
  }
  
  mu     = int.mod$fitted
  vari   = family$variance(mu)
  deriv  = 1/vari
  wt.mat = ob.wt*1/(deriv^2*vari)
  Xw.s   = t(as.vector(wt.mat) * t(t(X.des)))
  score  = t(as.vector(deriv) * t(Xw.s)) %*% (y - mu)
  score
}

singleLasso = function(y, X, lamb, ob.wt = rep(1, length(y)), alpha = 1, b.init = NA, intercept = NA, family = "gaussian", tol = 1.e-9, gamma = 0, init.coef = NA, mu.min = 1.e-16, mu.max = 1/mu.min, area.int = FALSE, interactions, max.it = 25, standardise = TRUE)
{
  error.flag = FALSE
  
  eta.min = family$linkfun(mu.min)
  eta.max = family$linkfun(mu.max)
  
  if (gamma == 0)
  {
    adapt.weights = rep(1, dim(X)[2])
  }
  if (gamma != 0)
  {
    adapt.weights = 1/abs(init.coef)^gamma
  }
  
  b.glm   = rep(NA, dim(X)[2])
  b.lasso = b.glm
  if (length(lamb) == 1)
  {
    lambda.start    = rep(lamb, dim(X)[2])	
    lambda.start[1] = 0
    lambda = as.array(lambda.start)
    lambda = lambda * abs(adapt.weights[1:length(lambda)])
  }
  
  if (length(lamb) > 1)
  {
    lambda.start = rep(0,dim(X)[2])
    lambda = as.array(lambda.start)
    lambda = lambda * abs(adapt.weights[1:length(lambda)])
  }
  
  if (area.int == TRUE)
  {
    lambda = c(lambda, 0)
    X      = cbind(X, interactions)
    is.in  = rep(TRUE, dim(X)[2])
  }
  
  killed  = is.infinite(lambda)
  X       = X[,killed == FALSE]
  b.lasso = b.lasso[killed == FALSE]
  lambda  = lambda[killed == FALSE]
  
  if (length(lamb) == 1 & lamb[1] == 0)
  {
    b.init    = NA
    intercept = NA
  }
  
  if (any(is.na(b.init)) & is.na(intercept))
  {
    mu.old  = rep(mean(y), dim(X)[1])
    eta.old = etaFromMu(mu.old, family = family, mu.min, mu.max, eta.min, eta.max)
    l.old   = calculateLikelihood(X, family, ob.wt, mu = mu.old, y, alpha, lambda, beta = rep(0, dim(X)[2]), penalty = FALSE)
    b.old   = rep(1, dim(X)[2])
    diff    = 10
    
    while (diff > tol)
    {
      update   = irlsUpdate(y, X, ob.wt, is.in, signs, eta = eta.old, mu = mu.old, alpha = alpha, lambda = lambda, beta.old = b.old, penalty = FALSE, family = family, mu.min, mu.max, eta.min, eta.max, tol)
      l.new    = calculateLikelihood(X, family, ob.wt, mu = update$mu, y, alpha, lambda, beta = update$beta, penalty = FALSE)
      diff     = abs(l.new - l.old)
      mu.old   = update$mu
      eta.old  = update$eta
      b.old    = update$beta
      l.old    = l.new
    }
    
    b.glm   = update$beta
    b.lasso = b.glm
  }
  
  if (any(is.na(b.init)) == FALSE)
  {
    b.lasso = b.init[killed == FALSE]
  }
  
  if (is.na(intercept) == FALSE)
  {
    if (family$family == "poisson")
    {
      b.lasso = c(log(mean(y)), rep(0, (dim(X)[2] - 1 - area.int)), rep(1, area.int))
    }
    
    if (family$family == "binomial")
    {
      b.lasso = c(log(mean(y)/(1 - mean(y))), rep(0, (dim(X)[2] - 1 - area.int)), rep(1, area.int))
    }
    
    if (family$family == "gaussian")
    {
      b.lasso = c(mean(y), rep(0, (dim(X)[2] - 1 - area.int)), rep(1, area.int))
    }
  }
  
  is.in        = abs(b.lasso) > 1.e-7
  sign.change  = rep(-1, dim(X)[2])
  is.different = is.in
  
  diff   = 10
  num.it = 0
  
  signs = sign(b.lasso)
  
  betas   = c(b.lasso)
  scores  = c()
  viols   = c()
  likes   = c()
  actions = c()
  
  eta.old = X %*% b.lasso
  eta.old[eta.old < eta.min] = eta.min
  eta.old[eta.old > eta.max] = eta.max
  mu.old  = muFromEta(eta.old, family = family, mu.min, mu.max)
  
  like    = calculateLikelihood(X, family, ob.wt, mu = mu.old, y, alpha, lambda, beta = b.lasso, penalty = TRUE)
  likes   = c(likes, like)
  
  if (length(lamb) == 1 & lamb[1] == 0)
  {	
    diff   = 0
    num.it = max.it + 1
  }
  
  while(diff > tol & num.it < max.it)
  {
    mult     = 0
    update   = irlsUpdate(y, X, ob.wt, is.in, signs, eta = eta.old, mu = mu.old, alpha = alpha, lambda = lambda, beta.old = b.lasso, penalty = TRUE, family = family, mu.min, mu.max, eta.min, eta.max, tol)
    
    
    if (update$error == "Singular matrix")
    {
      error.flag = TRUE
      break
    }
    rawdiff = b.lasso[is.in] - update$beta
    
    sign.change                   = signs[is.in]*sign(update$beta)
    sign.change[1]                = 1
    sign.change[sign.change == 0] = 1
    score.beta                    = rep(0, length(b.lasso))
    score.beta[is.in]             = update$beta
    
    if (any(sign.change != 1) == TRUE)
    {
      delta                       = update$beta - b.lasso[is.in]
      prop                        = min(abs(b.lasso[is.in]/delta)[sign.change != 1])
      b.lasso[is.in]              = b.lasso[is.in] + prop*delta
      is.in[abs(b.lasso) < 1.e-7] = FALSE
      b.lasso[is.in == FALSE]     = 0
      signs                       = sign(b.lasso)
      mult                        = 1
      score.beta                  = b.lasso
      
      update$eta                        = as.matrix(X[,is.in]) %*% b.lasso[is.in]
      update$eta[update$eta < eta.min]  = eta.min
      update$eta[update$eta > eta.max]  = eta.max
      update$mu                         = muFromEta(update$eta, family = family, mu.min, mu.max)
    }
    
    score                   = t(as.vector(update$deriv) * t(update$s.wt)) %*% (y - update$mu)
    score.lamb              = alpha*lambda + (1 - alpha) * lambda * as.vector(score.beta)
    score.lamb[lambda == 0] = 100000
    viol                    = as.vector(abs(score))/as.vector(score.lamb)
    bigviol                 = max(viol)
    
    if (any(sign.change != 1) != TRUE)
    {
      b.lasso[is.in] = update$beta
      signs          = sign(b.lasso)
      if (bigviol > 1 + 1.e-6 & is.in[viol == bigviol] == FALSE)
      {
        is.in[viol == bigviol][1]        = TRUE
        is.different[viol == bigviol][1] = TRUE
        signs[viol == bigviol][1]        = sign(score[viol == bigviol])
        mult                             = 1
      }
    }
    
    like    = calculateLikelihood(X, family, ob.wt, mu = update$mu, y, alpha, lambda, beta = b.lasso, penalty = TRUE)
    mu.old  = update$mu
    eta.old = update$eta
    
    for (act in 1:length(sign(b.lasso)))
    {
      if (sign(b.lasso)[act] == 0 & sign(as.matrix(betas)[act,dim(as.matrix(betas))[2]]) != 0)
      {
        actions = rbind(actions, paste("Step ", dim(as.matrix(betas))[2] + 1, ": Delete variable ", act, sep = ""))
      }
      if (sign(b.lasso)[act] != 0 & sign(as.matrix(betas)[act,dim(as.matrix(betas))[2]]) == 0)
      {
        actions = rbind(actions, paste("Step ", dim(as.matrix(betas))[2] + 1, ": Add variable ", act, sep = ""))
      }
    }
    
    betas        = cbind(betas, b.lasso)
    scores       = cbind(scores, score)
    viols        = cbind(viols, viol)
    likes        = cbind(likes, like)
    diff         = mult + abs(likes[length(likes)] - likes[length(likes) - 1])
    num.it       = num.it + 1
  }
  
  if (error.flag == TRUE)
  {
    return(list(b = NA, mu = NA, e.df = NA, deviance = NA, likelihood = NA, GCV = NA, AIC = NA, BIC = NA, HQC = NA, AICc = NA, ls = NA, pen.likelihood = NA, bs = NA, s = NA, v = NA, actions = NA, flag = "Singular matrix"))
    cat(paste("Singular matrix error. No model fit.", "\n"))
    flush.console()
    stop
  }
  
  if (error.flag != TRUE)
  {
    eta          = as.matrix(X[,is.in]) %*% as.matrix(b.lasso[is.in])
    eta[eta < eta.min] = eta.min
    eta[eta > eta.max] = eta.max
    mu           = muFromEta(eta, family = family, mu.min, mu.max)
    if (area.int != FALSE)
    {
      is.in[dim(X)[2]] = FALSE
    }
    wt           = irlsUpdate(y, X, ob.wt, is.in, signs, eta, mu, alpha, lambda, beta.old = b.lasso, penalty = TRUE, family = family, mu.min, mu.max, eta.min, eta.max, tol)$wt
    xwx          = wt %*% X[,is.in]
    bpinv        = rep(0, length(b.lasso))
    bpinv[is.in] = 1/abs(b.lasso[is.in])
    ginv         = diag(bpinv)	
    n.param      = length(b.lasso[is.in]) - 1 + area.int
    
    if (sum(is.in) > 1)
    {
      eff.df = sum(diag(solve(xwx + diag(as.vector(lambda[is.in])) %*% ginv[is.in, is.in]) %*% xwx)) + area.int
    }
    if (sum(is.in) == 1)
    {
      eff.df = sum(diag(solve(xwx + diag(as.matrix(lambda[is.in])) %*% ginv[is.in, is.in]) %*% xwx)) + area.int
    }	
    
    n.pres = sum(y > 0)
    
    unp.likelihood = calculateLikelihood(X, family, ob.wt, mu, y, alpha, lambda, beta = b.lasso, penalty = FALSE)
    pen.likelihood = calculateLikelihood(X, family, ob.wt, mu, y, alpha, lambda, beta = b.lasso, penalty = TRUE)
    
    if (family$family == "poisson")
    {
      dev            = -2*(sum(ob.wt*(y*log(mu) - mu)) + sum(y > 0)) 
    }
    
    if (family$family == "binomial")
    {
      dev            = -2*unp.likelihood
    }
    
    if (family$family == "gaussian")
    {
      dev            = -2*(unp.likelihood - calculateLikelihood(X, family, ob.wt[y > 0], y[y > 0], y, alpha, lambda, beta = b.lasso, penalty = FALSE))
    }
    
    if (eff.df < n.pres)
    {
      gcv     = dev/(n.pres*(1 - (eff.df)/n.pres)^2)
    }
    if (eff.df >= n.pres)
    {
      gcv = NA
    }
    
    aic = -2*unp.likelihood + 2*(n.param + 1)
    bic = -2*unp.likelihood + log(n.pres)*(n.param + 1)
    hqc = -2*unp.likelihood + 2*(n.param + 1)*log(log(n.pres))
    aicc = aic + 2*(n.param + 1)*(n.param + 2)/(n.pres - n.param - 2)
    
    b.out = rep(NA, length(killed))
    b.out[which(killed == FALSE)] = b.lasso
    b.out[which(killed == TRUE)] = 0
    
    return(list(b = b.out, mu = mu, e.df = eff.df, deviance = dev, likelihood = unp.likelihood, GCV = gcv, AIC = aic, BIC = bic, HQC = hqc, AICc = aicc, ls = likes, pen.likelihood = likes[length(likes)], bs = betas, s = scores, v = viols, actions = actions))
    
  }
}

ppmSpatstat = function(fit)
{
  is.ai  = is.numeric(fit$pt.interactions)
  pres.x = fit$x[fit$pres > 0]
  pres.y = fit$y[fit$pres > 0]
  quad.x = fit$x[fit$pres == 0]
  quad.y = fit$y[fit$pres == 0]
  
  ux = unique(quad.x)
  uy = unique(quad.y)
  ux = sort(ux)
  uy = sort(uy)
  nx = length(ux)
  ny = length(uy)
  
  quad.mask = matrix(NA, ny, nx, dimnames = list(uy, ux))
  col.ref   = match(quad.x, ux)
  row.ref   = match(quad.y, uy)
  
  all.vec = rep(NA, max(row.ref)*max(col.ref))
  vec.ref = (col.ref - 1)*max(row.ref) + row.ref
  all.vec[vec.ref] = 1
  num.vec = all.vec
  num.vec[is.na(all.vec)] = 0
  
  quad.mask = matrix(all.vec, max(row.ref), max(col.ref), dimnames = list(uy, ux))
  quad.im   = im(quad.mask, xcol = ux, yrow = uy)
  quad.win  = as.owin(quad.im)
  
  pres.dat = ppp(pres.x, pres.y, window = quad.win, check = FALSE)
  quad.dat = ppp(quad.x, quad.y, window = quad.win, check = FALSE)
  Q        = quad(data = pres.dat, dummy = quad.dat, w = fit$wt, param = list(weight = list(method = "grid", ntile = c(length(ux), length(uy)), npix = NULL, areas = num.vec)))
  
  num.var = dim(fit$data)[2] - 1 - is.ai
  
  cov.list = vector('list', num.var)
  
  for (var in 1:num.var)
  {
    x.dat = fit$data[fit$pres == 0, (var + 1)]
    v.vec = rep(NA, max(row.ref)*max(col.ref))
    v.vec[vec.ref] = x.dat
    x.mat = matrix(v.vec, max(row.ref), max(col.ref), dimnames = list(uy, ux))
    assign(paste("var.im.", var, sep = ""), im(x.mat, xcol = ux, yrow = uy))
    cov.list[[var]] = im(x.mat, xcol = ux, yrow = uy)
  }
  
  names(cov.list) = names(fit$beta)[-1]
  trend = paste("~", names(cov.list)[1], sep = "")
  for (var in 2:num.var)
  {
    trend = paste(trend, " + ", names(cov.list)[var], sep = "")
  }
  
  glmfit.form = as.formula(paste(".mpl.Y", trend))
  if (is.ai)
  {
    glmfit.form = as.formula(paste(".mpl.Y", trend, " + Interaction"))
  }
  
  trend  = as.formula(trend)
  call.1 = quote(ppm(Q = Q, trend = trend, covariates = cov.list, interaction = Poisson(), correction = "none"))
  if (is.ai)
  {
    call.1 = bquote(ppm(Q = Q, trend = trend, covariates = cov.list, interaction = AreaInter(.(fit$r)), correction = "none"))
  }
  
  class(fit)       = "ppm"
  fit$fitter       = "glm"
  fit$coef         = fit$beta
  fit$method       = "mpl"
  fit$projected    = FALSE
  fit$trend        = trend
  fit$interaction  = NULL
  if (is.ai)
  {
    fit$interaction = AreaInter(fit$r)
  }
  fit.int          = list(name = "Poisson process", creator = "Poisson", family = NULL, pot = NULL, par = NULL, parnames = NULL)
  class(fit.int)   = "interact"
  fit$fitin        = list(interaction = Poisson(), coefs = fit$beta, Vnames = character(0), IsOffset = logical(0))
  if (is.ai)
  {
    fit$fitin        = list(interaction = AreaInter(fit$r), coefs = fit$beta, Vnames = "Interaction", IsOffset = FALSE)
  }
  class(fit$fitin) = c("fii", "list")
  fit$Q            = Q
  fit$maxlogpl     = fit$likelihood
  fit$covariates   = cov.list
  fit$covfunargs   = list()
  fit$correction   = "none"
  fit$rbord        = 0
  
  glmdata          = data.frame(fit$wt, fit$pres/fit$wt, fit$data[,-1], TRUE)
  if (is.ai == FALSE)
  {
    names(glmdata)   = c(".mpl.W", ".mpl.Y", names(cov.list), ".mpl.SUBSET")
  }
  if (is.ai)
  {
    names(glmdata)   = c(".mpl.W", ".mpl.Y", names(cov.list), "Interaction", ".mpl.SUBSET")
  }
  fit$version     = list(major = 1, minor = 31, release = 0)
  fit$problems    = list()
  fit$call        = call.1
  fit$callstring  = "character"
  fit$callframe   = environment()
  
  terms.int                 = terms(glmfit.form)
  terms.ppm                 = terms(trend)
  fit$terms                 = terms.ppm
  fit$internal$glmfit$terms = terms.int
  
  #fullform = as.formula(paste("fit$pres/fit$wt", paste(as.character(fit$formula), collapse = '')))
  fullform = as.formula(paste("fit$pres/fit$wt", paste(as.character(trend), collapse = '')))
  
  glm.1 = glm(fullform, data = data.frame(fit$data), weights = fit$wt, family = poisson())
  #glm.1 = glm(fit$pres/fit$wt ~ as.matrix(fit$data[,-1]), weights = fit$wt, family = poisson())
  glm.1$coefficients      = fit$beta
  glm.1$fitted.values     = fit$mu
  glm.1$residuals         = (fit$pres/fit$wt - fit$mu)/fit$mu
  glm.1$linear.predictors = log(fit$mu)
  
  fam   = poisson()
  vari  = fam$variance(glm.1$fitted.values)
  deriv = 1/vari
  wt    = fit$wt
  weii  = wt*1/(deriv^2*vari)
  w.1   = weii
  Xw    = t(as.vector(sqrt(weii)) * t(t(fit$data)))
  q.1   = qr(t(Xw))
  
  glm.1$qr      = q.1
  glm.1$weights = w.1
  glm.1$family  = quasi(log)
  glm.1$data    = glmdata
  glm.1$formula = glmfit.form
  
  fit$internal  = list(glmfit = glm.1, glmdata = glmdata)
  if (is.ai)
  {
    fit$internal  = list(glmfit = glm.1, glmdata = glmdata, Vnames = "Interaction", IsOffset = FALSE)
  }
  
  fit
}

weights = function(sp.xy, quad.xy, coord = c("X", "Y"))
{
  sp.col   = c(which(names(sp.xy) == coord[1]), which(names(sp.xy) == coord[2]))
  quad.col = c(which(names(quad.xy) == coord[1]), which(names(quad.xy) == coord[2]))
  
  X.inc   = sort(unique(quad.xy[,quad.col[1]]))[2] - sort(unique(quad.xy[,quad.col[1]]))[1]
  Y.inc   = sort(unique(quad.xy[,quad.col[2]]))[2] - sort(unique(quad.xy[,quad.col[2]]))[1]
  quad.0X = min(quad.xy[,quad.col[1]]) - floor(min(quad.xy[,quad.col[1]])/X.inc)*X.inc
  quad.0Y = min(quad.xy[,quad.col[2]]) - floor(min(quad.xy[,quad.col[2]])/Y.inc)*Y.inc
  
  X = c(sp.xy[,quad.col[1]], quad.xy[,quad.col[1]])
  Y = c(sp.xy[,quad.col[2]], quad.xy[,quad.col[2]])
  
  round.X     = round((X - quad.0X)/X.inc)*X.inc
  round.Y     = round((Y - quad.0Y)/Y.inc)*Y.inc
  round.id    = paste(round.X, round.Y)
  round.table = table(round.id)
  wt          = X.inc*Y.inc/as.numeric(round.table[match(round.id, names(round.table))])
  
  wt
}

print.ppmlasso = function(x, ..., output = c("all", "path", "model", "interaction"))
{
  if ("all" %in% output)
  {
    if (x$family == "poisson")
    {
      output = c("path", "model")
    }
    if (x$family == "area.inter")
    {
      output = c("path", "model", "interaction")
    }	
  }
  model.grammar = if (length(x$lambdas) > 1)
    "models"
  else "model"
  family.name   = if (x$family == "poisson")
    paste("Poisson point process", model.grammar)
  else paste("Area-interaction", model.grammar, "(r =", x$r, ")")
  single.model  = if (x$family == "poisson")
    "Poisson point process model"
  else paste("Area-interaction model (r =", x$r, ")")
  pen.type = "LASSO"
  if (x$gamma > 0)
  {
    pen.type = paste("Adaptive LASSO (gamma =", x$gamma, ")")
  }
  if (x$alpha < 1)
  {
    pen.type = paste("Elastic Net (alpha =", x$alpha, ")")
  }
  if ("path" %in% output)
  {
    cat(paste("Regularisation path of", length(x$lambdas), family.name), "\n")
    cat(paste(pen.type, "Penalties:\n"))
    print(x$lambdas)
    cat("\nCoefficients:\n")
    print(x$betas)
    cat("\nLikelihoods:\n")
    print(x$likelihoods)
    cat("\nPenalised Likelihoods:\n")
    print(x$pen.likelihoods)
  }
  if ("model" %in% output)
  {
    cat(paste("Optimal", single.model, "chosen by", x$criterion, ":\n"))
    cat(paste(pen.type, "Penalty:\n"))
    print(x$lambda)
    cat("\nCoefficients:\n")
    print(x$beta)
    cat("\nLikelihood:\n")
    print(x$likelihood)
  }
  if ("interaction" %in% output)
  {
    cat(paste("Point interactions of radius", x$r, ":\n"))
    print(x$pt.interactions)
  }
}

setClass("ppmlasso")
setMethod("print", "ppmlasso", print.ppmlasso)

envelope.ppmlasso = function(Y, fun = Kest, ...)
{
  fit.ss = ppmSpatstat(Y)
  envelope(Y = fit.ss, fun = fun, ...)
}

setMethod("envelope", "ppmlasso", envelope.ppmlasso)

diagnose = function(object, ...)
  UseMethod("diagnose")

diagnose.ppmlasso = function(object, ...)
{
  fit.ss = ppmSpatstat(object)
  diagnose.ppm(fit.ss, ...)
}

setClass("ppm")
setGeneric("diagnose")
setMethod("diagnose", "ppmlasso", diagnose.ppmlasso)
setMethod("diagnose", "ppm", diagnose.ppm)

predict.ppmlasso = function(object, ..., newdata, interactions = NA)
{
  if (any(lapply(newdata, class) == "factor"))
  {
    unpacknewdata = catConvert(newdata)
    newdata       = unpacknewdata$X
    cat.names     = setdiff(unique(unpacknewdata$cat.names), NA)
    use.form = as.character(object$formula)[2]
    for (i in 1:length(cat.names))
    {
      use.form = gsub(cat.names[i], paste(names(newdata)[which(unpacknewdata$cat.names == cat.names[i])], collapse = " + "), use.form)
    }
    object$formula = as.formula(paste("~", use.form))
  }
  var.0 = which(apply(newdata, 2, var) == 0)
  #if (length(var.0) > 0)
  #{
  #newdata[,var.0] = newdata[,var.0] + rnorm(dim(newdata)[1], 0, 1.e-8)
  #}
  mf    = model.frame(object$formula, data = newdata)
  mt    = attr(mf, "terms")
  X.des = if (!is.empty.model(mt)) 
    model.matrix(mt, mf)
  else matrix(, length(object$mu), 0L)
  X.var = X.des[,-1]
  if (is.null(object$s.means) == FALSE)
  {
    X.var = scale(X.var, center = object$s.means, scale = object$s.sds)
    X.des = cbind(1, X.var)
  }
  if (object$family == "area.inter")
  {
    if (is.na(interactions) == TRUE)
    {
      if (is.null(object$s.means) == FALSE)
      {
        X.des = cbind(X.des, min(scale(object$pt.interactions)))
      }
      if (is.null(object$s.means) == TRUE)
      {
        X.des = cbind(X.des, 0)
      }
    }
    if (is.na(interactions) == FALSE)
    {
      if (is.null(object$s.means) == FALSE)
      {
        X.des = cbind(X.des, scale(interactions, center = mean(object$pt.interactions), scale = sd(object$pt.interactions)))
      }
      if (is.null(object$s.means) == TRUE)
      {
        X.des = cbind(X.des, interactions)
      }
    }
  }
  pred.int = exp(as.matrix(X.des) %*% object$beta)
  return(pred.int)
}

setMethod("predict", "ppmlasso", predict.ppmlasso)

blocks = function(n.blocks, block.scale, dat, seed = 1)
{
  if (is(dat, "list"))
  {
    dat = dat$X
  }
  cell.group = rep(0, length(dat$X))
  n.groups   = ceiling(c(max((dat$X - min(dat$X))/block.scale), max((dat$Y - min(dat$Y))/block.scale)))
  
  xq = block.scale*(1:n.groups[1]) + min(dat$X)
  for (i.group in 1:n.groups[1])
  {
    cell.group = cell.group + as.numeric(dat$X > xq[i.group])
  }
  
  yq = block.scale*(1:n.groups[2]) + min(dat$Y)
  for (i.group in 1:n.groups[2])
  {
    cell.group = cell.group + n.groups[1] * as.numeric(dat$Y > yq[i.group])
  }
  
  block.group = factor(cell.group)
  
  set.seed(seed) #to ensure that you get the same sample each time - useful for back-tracking
  
  levs       = rep(1:n.blocks, length = length(levels(block.group)))
  lev.sample = sample(levs)
  
  levels(block.group) = lev.sample
  block.group
}

calculateUnpenalisedLikelihood = function(ob.wt, y, mu)
{
  like = sum(ob.wt*(y*log(mu) - mu)) - sum(log(1:sum(y > 0)))
  like
}
