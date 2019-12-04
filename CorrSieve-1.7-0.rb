#-----------------------------------------------------------------------------------------------
# CorrSieve 1.7-0: software to summarise and evaluate output from STRUCTURE
# Noted code improvements by Michael G. Campana (2019) are in the public domain
# as a US government work.
#
# Original source code is from:
# CorrSieve 1.6-5: software to summarise and evaluate output from STRUCTURE
# Copyright (C) 2010-2011 Michael G. Campana
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------------------------

class Array
#-----------------------------------------------------------------------------------------------
	# Method each_permutation removed in CorrSieve 1.7.0 update
	def inexact
		a = []
		b = []
		q = self.size - 1
		for count in 0 .. q
			ran = rand(q)
			while b.include?(count)
				ran += 1
				ran = 0 if ran > q
			end
			b.push(count)
			a.push(self[ran])  
		end
		return a
	end
end
#------------------------------------------------
# Universal Variables
#------------------------------------------------
$k = 1 # maximum K value
$numK = [] # iterations of K
$save = "" # save file path
$corr = false # calculate matrix correlations
$r2 = 0.7 # minimum r for filter
$usep = false # determines whether to use permutation test for p
$p = 0.05 # minimum p for filter
$maxcorr = false
$input_arrays = [] # input data array
$print = false # option to print correlations matrices
$inexact = false # use inexact p
$iterat = 1 # number of iterations to estimate p
$deltaK = false # whether to calculate deltaK
$dKarray = [] # array to store LnPD data
$Fst = false # calculate mean Fsts
$Fstarray = [] # array to store Fst data
$Fstopt = 1 #method for Fst statistics optimisation 
#-------------------------------------------------
# Additional Methods
#-------------------------------------------------
def mean(val = [])
	mean = 0.0
	for i in 0.. val.size-1
		mean += val[i]
	end
	mean /= val.size
	return mean
end
#-------------------------------------------------
def stdev(val = [])
	me = mean(val)
	st = 0.0
	for stval in val
		add = (stval - me) * (stval - me)
		st += add
	end
	de = (val.size - 1).to_f
	st /= de
	st2 = Math.sqrt(st)
	return st2
end
#-------------------------------------------------
def calc_dK
	dk = "K\tMean lnP(D)\tL'(K)\tL''(K)\tdelta K\n"
	means = []
	stdevs = []
	first = []
	second = []
	dKs = []
	for i in 1 .. $k
		if $numK[i-1] > 1
			me = mean($dKarray[i-1])
			means.push(me)
			stdevs.push(stdev($dKarray[i-1]))
		else
			means.push("NA")
			stdevs.push("NA")
		end
	end
	for i in 0 .. $k - 1
		if (means[i] != "NA" and means[i-1] != "NA" and i != 0)
			val = means[i]-means[i-1]
			first.push(val)
		else
			first.push("NA")
		end
	end
	for i in 0 .. $k - 1
		if (first[i+1] != "NA" and first[i] != "NA" and stdevs[i] != "NA" and i != $k - 1)
			val = first[i+1] - first[i]
			second.push(val.abs)
			dKs.push(val.abs/stdevs[i])
		else
			second.push("NA")
			dKs.push("NA")
		end
	end
	for i in 1 .. $k
		add = i.to_s + "\t" + means[i-1].to_s + "\t" + first[i-1].to_s + "\t" + second[i-1].to_s + "\t" + dKs[i-1].to_s + "\n"
		dk+= add
	end
	return dk
end
#-------------------------------------------------
def deltaF(df = [])
	first = []
	dFs = []
	means = []
	stdevs = []
	for i in 1 .. $k
		if $numK[i-1] > 1
			means.push(mean(df[i-1]))
			stdevs.push(stdev(df[i-1]))
		else
			means.push("NA")
			stdevs.push("NA")
		end
	end
	for i in 0 .. $k - 1
		if (means[i] != "NA" and means[i-1] != "NA" and i != 0)
			first.push(means[i] - means[i -1])
		else
			first.push("NA")
		end
	end
	for i in 0 .. $k - 1
		if (first[i+1] != "NA" and first[i] != "NA" and stdevs[i] != "NA" and i != $k - 1)
			dFs.push((first[i+1] - first[i]).abs/stdevs[i])
		else
			dFs.push("NA")
		end
	end
	return dFs
end
#-------------------------------------------------
def pvalue(roar, alpha = [], beta = [])
	greater = 0.0
	total = 0.0
	perm = alpha.permutation.to_a # Improvement of permutation code in public domain [line 164]
	for i in 0 .. perm.size-1
		greater += 1.0 if corr(perm[i], beta).abs >= roar
		total += 1.0
	end
	p = greater/total
	return p
end
#-------------------------------------------------
def corr(alpha = [], beta =[])
	meanx = mean(alpha)
	meany = mean(beta)
	numer = 0.0
	denomx = 0.0
	denomy = 0.0
	for i in 0..alpha.size-1
		numer += ((alpha[i]-meanx) * (beta[i]-meany))
		denomx += ((alpha[i]-meanx)*(alpha[i]-meanx))
		denomy += ((beta[i]-meany) * (beta[i]-meany))
	end
	r = numer/(Math.sqrt(denomx)*Math.sqrt(denomy))
	return r
end 
#-------------------------------------------------
def import_data(file)
	kbreak = false
	start = false
	data = []
	count = 0
	fstopt = []
	File.open(file, 'r') do |imp|
		while line = imp.gets
			break if kbreak
			if line.include?("populations assumed")
				line2 = line.delete " populations assumed"
				@kvalue = line2.to_i
				$numK[@kvalue - 1] += 1
				kbreak = true
			elsif line.include?("Population number assumed=")
				line2 = line.delete "Population number assumed="
				@kvalue = line2.to_i
				$numK[@kvalue - 1] += 1
				kbreak = true
			end
		end
	end
	for i in 1 .. @kvalue
		data.push([])
	end
	File.open(file, 'r') do |imp|
		while line = imp.gets
			break if start and line == "\n"
			if $Fst
				if line.include?("Mean value of Fst_")
					if $Fstopt == 2
						fstopt.push(line[29+(count+1).to_s.length,6].to_f)
					else
						$Fstarray[@kvalue - 1][count].push(line[29+(count+1).to_s.length,6].to_f)						
					end
					count += 1
				end
			end
			if $deltaK
				if line.include?("Estimated Ln Prob of Data")
					line2 = line.delete "Estimated Ln Prob of Data   = "
					$dKarray[@kvalue - 1].push line2.to_f
				elsif line.include?("Posterior Mean = ")
					line2 = line.delete "Posterior Mean = "
					$dKarray[@kvalue - 1].push line2.to_f
				end
			end
			if (start and ($corr or $Fstopt == 3))
				while line.include?(":")
					line[0,1] = ""
				end 
				line.gsub!("\t"," ") if line.include?("\t")
				line[0,2] = ""
				for i in 0 .. @kvalue - 1
					val = line[i * 6,5]
					data[i].push val.to_f
				end
			end
			start = true if line.include?("Label (%Miss)")
			start = true if line.include?("Label\t(Miss)")
		end
	end
  	if $Fstopt == 2
		fstopt.sort!
		for x in 0 .. count - 1
			$Fstarray[@kvalue - 1][x].push(fstopt[x])
		end
	end
	$input_arrays[@kvalue - 1].push(data)
end
#-------------------------------------------------
def valid(val = "")
	while (val.upcase != "Y" and val.upcase != "YES" and val.upcase != "N" and val.upcase != "NO")
		puts "Invalid response. Please re-enter (Y/N)."
		val = gets.chomp
	end
	return val
end
#-------------------------------------------------
# Processing Block: Calculations begin here. 
#-------------------------------------------------
def process_data # Conversion to method in CorrSieve 1.7-0 update (Lines 269-270 is public domain)
	puts "Please enter the path to the file folder."
	@folder = gets.chomp
	while !FileTest.directory?(@folder)
		print "Folder not found. Please re-enter the path carefully.\n"
		@folder = gets.chomp
	end
	puts "Please enter path for the saved files.\nPressing ENTER will keep the current directory."
	$save = gets.chomp
	$save = File.expand_path("") if $save == ""
	Dir.mkdir($save) if !FileTest.directory?($save)
	@run = ""
	while @run == ""
		puts "Please enter the run name."
		@run = gets.chomp
	end
	$save  += "/"
	$save += @run
	puts "Calculate Q matrix correlations? (Y/N)"
	tmp = gets.chomp
	tmp = valid(tmp)
	if tmp.upcase == "Y" or tmp.upcase == "YES"
		$corr = true
	end
	if $corr
		puts "Please enter minimum r value.\n"
		$r2 = gets.chomp.to_f
		while $r2.abs > 1.0
			print "Invalid r value. Please re-enter.\n"
			$r2 = gets.chomp.to_f
		end
		puts "Use p-value cut-off? (Y/N)"
		tmp = gets.chomp
		tmp = valid(tmp)
		if tmp.upcase == "Y" or tmp.upcase == "YES"
			$usep = true
			puts "Estimate p? (Y/N). If no, CorrSieve will calculate an\nexact p. This will be very slow."
			tmp = gets.chomp
			tmp = valid(tmp)
			if tmp.upcase == "Y" or tmp.upcase == "YES"
				$inexact = true
				puts "Please enter the number of permutations to estimate p values."
				$iterat = gets.chomp.to_i
				while $iterat < 1
					print "Invalid number of permutations. Please re-enter.\n"
					$iterat = gets.chomp.to_i
				end
			end
			puts "Please enter maximum p value.\n"
			$p = gets.chomp.to_f
			while $p < 0.0 or $p > 1.0
				print "Invalid p value. Please re-enter.\n"
				$p = gets.chomp.to_f
			end
		end
		puts "Use average max correlation criterion? (Y/N). If no\nCorrSieve will use column and row criterion."
		tmp = gets.chomp
		tmp = valid(tmp)
		if tmp.upcase == "Y" or tmp.upcase == "YES"
			$maxcorr = true
		end
		puts "Output unfiltered correlation matrices? (Y/N)"
		tmp = gets.chomp
		tmp = valid(tmp)
		if tmp.upcase == "Y" or tmp.upcase == "YES"
			$print = true
		end
	end
	puts "Summarise Ln P(D) and calculate delta K? (Y/N)"
	tmp = gets.chomp 
	tmp = valid(tmp)
	if tmp.upcase == "Y" or tmp.upcase == "YES"
		$deltaK = true
	end
	puts "Calculate Fst statistics? (Y/N)"
	tmp = gets.chomp
	tmp = valid(tmp)
	if tmp.upcase == "Y" or tmp.upcase == "YES"
		$Fst = true
		puts "Select Fst statistics optimisation method:\n1. No optimisation\n2. Order Fst data by value\n3. Use matrix correlations"
		$Fstopt = gets.chomp.to_i
		while $Fstopt != 1 and $Fstopt !=2 and $Fstopt !=3
			puts "Invalid selection. Select Fst statistics optimisation method:\n1. No optimisation\n2. Order Fst data by value\n3. Use matrix correlations"
			$Fstopt = gets.chomp.to_i
		end
	end
  #-------------------------------------------------------
  # Data import block
  #-------------------------------------------------------
	puts "Processing...This may take a while. Tea break?"
	Dir.foreach(@folder + "/") do |imp|
		if imp[-2,2] == "_f"
			kbreak = false
			File.open(@folder + "/" + imp, 'r') do |kcheck|
				while line = kcheck.gets
					break if kbreak
					if line.include?("populations assumed")
						line2 = line.delete " populations assumed"
						$k = line2.to_i if $k < line2.to_i
						kbreak = true
					elsif line.include?("Population number assumed=")
						line2 = line.delete "Population number assumed="
						$k = line2.to_i if $k <line2.to_i
						kbreak = true
					end
				end
			end
		end
	end
	for i in 1 .. $k
		$numK.push(0)
		$input_arrays.push([])
		$dKarray.push [] if $deltaK
		if $Fst
			$Fstarray.push []
			for j in 1 .. $k - ($k - i)
				$Fstarray[i-1].push []
			end
		end
	end
	Dir.foreach(@folder + "/") do |imp|
		import_data(@folder + "/" + imp) if imp[-2,2] == "_f"
	end
	GC.start
  #-------------------------------------------------------
  # Data processing block
  # Vars:
  # i: determines the K value. This needs to iterate once for each K and will not permute
  # j: determines the run number of matrix 1. Needs to COMBINE
  # k: determines the col. of matrix 1. Needs to permute
  # l: determines the run number of matrix 2. Needs to COMBINE
  # m: determines the col. of matrix 2. Needs to permute
  #-----------------------------------------------------
	if $corr
		out = ""
		filter = ""
		if $maxcorr
			mcorr = 0
			mcarr = []
			mcfilt = "\n"
		else
			filtx = false
			filty = true
		end
		for i in 0 .. $k - 1
			if $numK[i] > 1
				add = "K = " + (i+1).to_s + "\n"
				out += add
				filter += add
				if $maxcorr
					filter += "Average Max Correlation = "
				else
					filter += " "
					for a in 2 .. $numK[i]
						add = "\t" + a.to_s
						filter += add
					end
					filter += "\n"
				end
				for j in 0 .. $numK[i] - 2
					la = j
					if !$maxcorr
						filter += (j+1).to_s
						if j > 0
							for f in 1 .. $k - ($k - j)
								filter += "\t"
							end
						end
					end
					for l in 1 .. $numK[i] - (j + 1)
						la += 1
						add = "Run " + (j+1).to_s + " versus Run " + (la+1).to_s + "\n"
						out += add
						mat = []
						matp = []
						for k in $input_arrays[i][j]
							x = []
							xp = []
							for m in $input_arrays[i][la]
								r = corr(k, m)
								x.push(r)
								out += r.to_s
								filtx = true if (r.abs >= $r2.abs and !$maxcorr)
								if $usep
									if $inexact
										greater = 0.0
										total = 0.0
										for u in 1 .. $iterat
											k1 = k.inexact
											m1 = m.inexact
											total += 1.0
											greater += 1.0 if corr(k1, m1).abs >= r
										end
										p = greater/total
										xp.push(p)
									else
										p = pvalue(r.abs,k,m)
										xp.push(p)
									end
									add = "(" + p.to_s + ")"
									out += add
									if $maxcorr
										mcorr = r.abs if (r.abs > mcorr and p <= $p)
									else
										filtx = false unless (r.abs >= $r2.abs and p <= $p)
									end
								elsif $maxcorr
									mcorr = r.abs if (r.abs > mcorr)
								end
								out += " "
								mat.push(x)
								matp.push(xp)
							end
							testy = []
							out += "\n"
							if $maxcorr
								mcarr.push mcorr
								mcorr = 0
							else
								for z in 0 .. mat.size - 1
									tempy = false
									for y in 0 .. mat[z].size - 1
										tempy = true if mat[z][y].abs >= $r2.abs
										if $usep
											tempy = false unless (mat[z][y].abs >= $r2.abs and matp[z][y].abs <= $p)
										end
									end
									testy.push(tempy)
								end
							end
						end
						if !$maxcorr
							filty = false if testy.include?(false)
							if filtx and filty
								filter += "\tY"
							else
								filter += "\tN"
							end
							filtx = false
							filty = true
						end
						out += "\n"
					end
					if !$maxcorr
						filter += "\n"
					end
				end
				out += "\n"
				if $maxcorr
					filter += mean(mcarr).to_s
					if mean(mcarr) >= $r2
						mcfilt += "\tY"
					else
						mcfilt += "\tN"
					end
				end
				mcarr = []
				filter += "\n\n"
			end
		end
		if $maxcorr
			filter += "Significant Correlations\nK = "
			for i in 1 .. $k
				if $numK[i-1] > 1
					filter += "\t"
					filter += i.to_s
				end
			end
			filter += mcfilt
		end
		if $print
			filename = $save + "-matrix_correlations.txt"
			File.open(filename, 'w') do |f2|
				f2.puts out
			end
		end
		filename = $save +"-filtered.txt"
		File.open(filename, 'w') do |f2|
			f2.puts filter
		end
	end
	if $deltaK
		filename = $save + "-deltaK.txt"
		File.open(filename, 'w') do |f2|
			f2.puts calc_dK
		end
	end
	if $Fst
		test = ""
		df = []
		dats = "\nK\tOverall Mean Fst\tOverall St. Dev.\tSt. Dev. of Means\tMean St. Dev.\tSt. Dev. of St. Devs.\tdelta Fst\n"
		datsarr = []
		for i in 1 .. $k
			tmp = []
			tmp2 = []
			tmp3 = []
			if $Fstopt == 3 and i > 1
				for j in 1 .. $numK[i-1] - 1
					xmat = []
					xtmp = []
					fstarr = []
					for k in $input_arrays[i-1][0]
						x = []
						for m in $input_arrays[i-1][j]
							r = corr(k, m)
							x.push(r)
						end
						xmat.push(x)
					end
					for f in 0 .. i - 1
						fstarr.push($Fstarray[i-1][f][j])
					end
					for f in 0 .. i -1
						rv = -1.0
						sp = 1
						for g in 0 .. i -1
							if xmat[f][g] > rv
								rv = xmat[f][g]
								sp = g
							end
						end
						xtmp.push(fstarr[sp])
						for h in xmat
							h[sp] = - 2
						end
					end
					for f in 0 .. i - 1
						$Fstarray[i-1][f][j] = xtmp[f]
					end
				end
			end
			if $numK[i-1] > 0
				add = "K = " + i.to_s + "\nCluster\tMean Fst\tSt. Dev.\n"
				test += add   
				for j in 1 .. $k - ($k - i)
					if $numK[i-1] > 1
						add = j.to_s + "\t" + mean($Fstarray[i-1][j-1]).to_s + "\t" + stdev($Fstarray[i-1][j-1]).to_s + "\n"
						tmp2.push(stdev($Fstarray[i-1][j-1]))
					else
						add = j.to_s + "\t" + $Fstarray[i-1][j-1].to_s + "\tNA\n"
						tmp2.push("NA")
					end
					for dat in $Fstarray[i-1][j-1]
						tmp3.push(dat)
					end
					tmp.push(mean($Fstarray[i-1][j-1]))
					test += add
				end
			else
				tmp.push("NA")
				tmp2.push("NA")
				tmp3.push("NA")
			end
			df.push(tmp3)
			if tmp.include?("NA")
				metmp = "NA"
				sdtmp = "NA"
			else
				metmp = mean(tmp).to_s
				if tmp.size > 1
					sdtmp = stdev(tmp).to_s
				else
					sdtmp = "NA"
				end
			end
			if tmp2.include?("NA")
				metmp2 = "NA"
				sdtmp2 = "NA"
			else
				metmp2 = mean(tmp2).to_s
				if tmp2.size > 1
					sdtmp2 = stdev(tmp2).to_s
				else
					sdtmp2 = "NA"
				end
			end
			if tmp3.include?("NA")
				metmp3 = "NA"
				sdtmp3 = "NA"
			else
				metmp3 = mean(tmp3).to_s
				if tmp3.size > 1
					sdtmp3 = stdev(tmp3).to_s
				else
					sdtmp3 = "NA"
				end
			end
			add = i.to_s + "\t" + metmp3 + "\t" + sdtmp3 + "\t" + sdtmp + "\t" + metmp2 + "\t" + sdtmp2 + "\t"
			datsarr.push(add)
		end        
		delFst = deltaF(df)
		for i in 1 .. datsarr.size
			add = datsarr[i-1] + delFst[i-1].to_s + "\n"
			dats += add
		end
		test += dats
		filename = $save +"-Fst.txt"
		File.open(filename, 'w') do |f2|
			f2.puts test
		end
	end
end
# Methods show_license [lines 671-695, 709-722], gnu_license [lines 724-737], and splash_screen [lines 741-763] are in the public domain. GNU General Public License header [lines 696-707] copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
#------------------------------------------------- 
# License Information
#-------------------------------------------------
def show_license
	if RUBY_PLATFORM =~ /win32/
		system('cls')
	else
		system("clear")
	end
	puts <<BASIC_LICENSE
CorrSieve 1.7-0.
Public domain updates by Michael G. Campana, 2019
Original CorrSieve 1.6-5 source code copyright (c) Michael G. Campana, 2010-2011
Code improvements by Michael G. Campana (2019) are US government works
and are therefore in the public domain. Affected lines of source code
are annotated that they are in the public domain.

Improvements include:
* Removal of the method 'each_permutation' from the Array class.
* Improvement of the 'pvalue' method.
* Conversion of the main processing section to the method 'process_data'.
* Revisions to the splash screen for updated licensing information.
* Addition of method 'splash_screen' to control revised program.
* Addition of methods 'show_license' and 'gnu_license' to detail licensing information.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Enter 'L' to see full GNU General Public License or 'R' to return to the splash screen.
BASIC_LICENSE
	choice = gets.chomp.upcase
	while choice != "R" && choice != "L"
		puts "\Enter 'L' to see full GNU General Public License or 'R' to return to the splash screen."
		choice = gets.chomp.upcase
	end
	case choice
	when 'L'
		gnu_license
	when 'R'
		splash_screen
	end
end
#-------------------------------------------------
def gnu_license
	if FileTest.exist?(File.absolute_path(__dir__) + "/LICENSE")
		File.open(File.absolute_path(__dir__) + "/LICENSE") do |f1|
			while line = f1.gets
				puts line 
			end
		end
	else
		puts "LICENSE file not found in CorrSieve executable's directory"
	end
	puts "\nPress Enter to return to the splash screen."
	cont = gets.chomp
	splash_screen
end
#-------------------------------------------------
# Options Block: Program begins here
#-------------------------------------------------
def splash_screen
	if RUBY_PLATFORM =~ /win32/
		system('cls')
	else
		system("clear")
	end
	print "Welcome to CorrSieve 1.7-0\n\nPublic domain updates by Michael G. Campana, 2019 (Smithsonian Conservation Biology Institute)\nOriginal CorrSieve source code copyright (c) Michael G. Campana, 2010-2011\n\nThis software is licensed under the GNU General Public License (version 3 or later)\nFor details, see the license provided with this code or at <http://www.gnu.org/licenses/>.\nThis program comes with absolutely no warranty.\n\n"
	puts "Enter 'C' to continue, 'L' to see license, or 'X' to exit"
	choice = gets.chomp.upcase
	while choice != "C" && choice != "L" && choice != "X"
		puts "Enter 'C' to continue, 'L' to see license details, or 'X' to exit"
		choice = gets.chomp.upcase
	end
	case choice
	when 'C'
		process_data
	when 'L'
		show_license
	when 'X'
		exit
	end
end
splash_screen