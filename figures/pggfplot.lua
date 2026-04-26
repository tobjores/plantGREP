-- function to generate a 'pggfplot' package error in LaTeX
local function pggf_error(message)
	tex.sprint('\\PackageError{pggfplots}{' .. message .. '}{}')
	error(message)
end


-- function to generate a 'pggfplot' package warning in LaTeX
	--[[ 
		if the optional `ifswitch` is given, the warning will only be issued if `\ifpggf@<ifswitch>` is true
	--]]
local function pggf_warning(message, ifswitch)
	if ifswitch then
		tex.sprint(4, '\\ifpggf@' .. ifswitch)
	end
	tex.sprint(4, '\\pggf@warning{' .. message .. '}')
	if ifswitch then
		tex.sprint('\\fi')
	end
end


-- function to test if a file exists
local function file_exists(file)
   local f = io.open(file, 'r')
   if f ~= nil then
		io.close(f)
		return true
	else
		return false
	end
end


-- function to remove the nth element from a string x
	--[[
		elements must be separated by tabs
	--]]
local function remove_nth_element(x, n)
	local new_string

	local i = 1

	for element in string.gmatch(x, '[^\t]*') do
		if i ~= n then
			if new_string then
				new_string = new_string .. '\t' .. element
			else
				new_string = element
			end
		end
		i = i + 1
	end

	return new_string
end


-- function to multiply the nth element from a string x with factor f and then add number a
	--[[
		elements must be separated by tabs
	--]]
local function modify_nth_element(x, n, f, a)
	local new_string

	local i = 1

	for element in string.gmatch(x, '[^\t]*') do
		if i == n then
			element = element * f + a
		end

		if new_string then
			new_string = new_string .. '\t' .. element
		else
			new_string = element
		end

		i = i + 1
	end

	return new_string
end


-- function to switch to elements n and m in a string x
	--[[
		elements must be separated by tabs
	--]]
local function switch_elements(x, n, m)
	local new_string
	local replaced

	local i = 1
	local element_n
	local element_m

	for element in string.gmatch(x, '[^\t]*') do
		if i == n then
			element_n = element
			element = element_m or '__placeholder_n__'
		elseif i == m then
			element_m = element
			element = element_n or '__placeholder_m__'
		end

		if new_string then
			new_string = new_string .. '\t' .. element
		else
			new_string = element
		end

		i = i + 1
	end

	new_string, replaced = string.gsub(new_string, '__placeholder_n__', element_m)

	if replaced == 0 then
		new_string  = string.gsub(new_string, '__placeholder_m__', element_n)
	end

	return new_string
end


-- function to get number of facet columns and rows
local function get_facet_ns(facet_table)
	local facets = {}

	for _, value in pairs(facet_table) do
		facets[value] = (facets[value] or 0) + 1
	end

	local dim1 = 0
	local dim2 = 0

	for key, value in pairs(facets) do
		dim1 = dim1 + 1
		dim2 = math.max(dim2, value)
	end

	return dim1, dim2
end


-- function to get facet info
local function get_facet_info(file)
	if not file_exists(file) then
		pggf_error('Could not read file \'' .. file .. '\'')
	end

	-- parse file
	local lines = io.lines(file)

	local header = lines()

	local columns = {}
	local facet_info = {}

	for colName in string.gmatch(header, '[^\t]*') do
		table.insert(columns, colName)
		table.insert(facet_info, {})
	end

	local col_to_id = {}

	for i, v in ipairs(columns) do
		col_to_id[v] = i
	end

	for line in lines do
		local i = 1
		for col in string.gmatch(line, '[^\t]*') do
			table.insert(facet_info[i], col)
			i = i + 1
		end
	end

	-- handle facet skipping
	local includefacets = token.get_macro('pggf@includefacets')
	local skipfacets = token.get_macro('pggf@skipfacets')

	if not string.find(includefacets, '^[0-9,%s]*$') then
		pggf_error('Invalid facet specification in `only facets = {' .. includefacets .. '}`. Can contain only comma-separated numbers')
	end
	if not string.find(skipfacets, '^[0-9,%s]*$') then
		pggf_error('Invalid facet specification in `skip facets = {' .. skipfacets .. '}`. Can contain only comma-separated numbers')
	end

	includefacets = string.gsub(includefacets, '%s+', '')
	skipfacets = string.gsub(skipfacets, '%s+', '')

	local facet_list = {}

	if includefacets == '' then
		for facet = 1, #facet_info[1] do
			if not string.find(',' .. skipfacets .. ',', ',' .. facet .. ',') then
				table.insert(facet_list, tostring(facet))
			end
		end
	else
		for facet in string.gmatch(includefacets, '[^,]+') do
			if not string.find(',' .. skipfacets .. ',', ',' .. facet .. ',') then
				table.insert(facet_list, facet)
			end
		end
	end

	if col_to_id['col'] and col_to_id['row'] then
		-- if plot contains rows and columns, skip only complete rows or columns
		local skipped_facets = {}
		local included_rows = {}
		local included_cols = {}
		local included_facets = {}

		for _, facet in ipairs(facet_list) do
			included_facets[facet] = true
		end

		for i = 1, #facet_info[1] do
			if included_facets[facet_info[col_to_id['facet_id']][i]] then
				included_cols[facet_info[col_to_id['col']][i]] = true
				included_rows[facet_info[col_to_id['row']][i]] = true
			else
				skipped_facets[facet_info[col_to_id['facet_id']][i]] = {}
				skipped_facets[facet_info[col_to_id['facet_id']][i]]['col'] = facet_info[col_to_id['col']][i]
				skipped_facets[facet_info[col_to_id['facet_id']][i]]['row'] = facet_info[col_to_id['row']][i]
			end
		end

		for f, cr in pairs(skipped_facets) do
			if included_cols[cr['col']] and included_rows[cr['row']] then
				table.insert(facet_list, f)
				pggf_warning('Cannot skip facet ' .. f .. ' because there are other facets in the same row and/or column')
			end
		end
	end

	table.sort(facet_list, function(a, b) return tonumber(a) < tonumber(b) end)

	FACETS = {}
	local p = 1

	for _, facet in ipairs(facet_list) do
		if not FACETS[facet] then
			FACETS[facet] = p
			p = p + 1
		end
	end

	local delete_rows = {}

	for i, v in ipairs(facet_info[col_to_id['facet_id']]) do
		-- delete skipped facets and rename others
		if FACETS[v] then
			facet_info[col_to_id['facet_id']][i] = FACETS[v]
		else
			table.insert(delete_rows, i)
		end
	end

	table.sort(delete_rows, function(a, b) return a > b end)

	for _, v in ipairs(delete_rows) do
		for k, _ in pairs(facet_info) do
			table.remove(facet_info[k], v)
		end
	end

	-- export number of columns/rows
	local n_cols
	local n_rows

	if col_to_id['col'] then
		n_cols, n_rows = get_facet_ns(facet_info[col_to_id['col']])
	elseif col_to_id['row'] then
		n_rows, n_cols = get_facet_ns(facet_info[col_to_id['row']])
	elseif col_to_id['wrap'] then
		n_cols = tonumber(token.get_macro('pggf@wrap@cols'))

		if n_cols < 1 then
			pggf_error('The number of columns has to be indicated with `wrap cols = <n>` for wrapped facets')
		end

		n_rows = math.ceil(#facet_info[1] / n_cols)

		token.set_macro('pggf@wrap@facets', #facet_info[1])
	else
		pggf_error('Could find neither \'col\', \'row\', nor \'wrap\' columns  in file \'' .. file .. '\'')
	end

	token.set_macro('pggf@facetcols', n_cols)
	token.set_macro('pggf@facetrows', n_rows)

	-- export name of column/row
	if col_to_id['col'] then
		for i, v in ipairs(facet_info[col_to_id['col']]) do
			tex.sprint('\\pgfplotsset{facet ' .. i .. '/.append style={/pggf/thiscolname={'.. v ..'}}}')
		end
	end

	if col_to_id['wrap'] then
		tex.sprint(4, '\\pggf@facet@wraptrue')

		for i, v in ipairs(facet_info[col_to_id['wrap']]) do
			tex.sprint('\\pgfplotsset{facet ' .. i .. '/.append style={/pggf/thiscolname={'.. v ..'}}}')
		end
	else
		tex.sprint(4, '\\pggf@facet@wrapfalse')
	end

	if col_to_id['row'] then
		for i, v in ipairs(facet_info[col_to_id['row']]) do
			tex.sprint('\\pgfplotsset{facet ' .. i .. '/.append style={/pggf/thisrowname={'.. v ..'}}}')
		end
	end

	-- export column/row titles
	if (n_cols > 1) or (col_to_id['wrap'] and (n_rows > 1)) then
		local col_titles = token.get_macro('pggf@coltitlesfromtable')

		if col_titles ~= 'false' then
			local title_col = col_to_id['col']

			if col_to_id['wrap'] then
				title_col = col_to_id['wrap']
			end

			local titles = facet_info[title_col][1]

			local n_titles = n_cols
			if col_to_id['wrap'] then
				n_titles = #facet_info[1]
			end

			for i = 2, n_titles do
				titles = titles .. ',' .. facet_info[title_col][i]
			end

			if col_titles == 'true' then
				token.set_macro('pggf@coltitles', titles)
			else
				token.set_macro(col_titles, titles, 'global')
			end
		end
	end

	if n_rows > 1 and not col_to_id['wrap'] then
		local row_titles = token.get_macro('pggf@rowtitlesfromtable')

		if row_titles ~= 'false' then
			local titles = facet_info[col_to_id.row][1]

			for i = 2 * n_cols, n_cols * n_rows, n_cols do
				titles = titles .. ',' .. facet_info[col_to_id.row][i]
			end

			if row_titles == 'true' then
				token.set_macro('pggf@rowtitles', titles)
			else
				token.set_macro(row_titles, titles, 'global')
			end
		end
	end
end


-- function to split a data file by facets
local function split_table_by_facets(file)
	if not file_exists(file) then
		pggf_error('Could not read file \'' .. file .. '\'')
	end

	local lines = io.lines(file)

	local header = lines()

	local tables = {}

	if string.find('\t' .. header .. '\t', '\tfacet_id\t') then
		local fid_col = 1
		for colName in string.gmatch(header, '[^\t]*') do
			if colName == 'facet_id' then
				break
			end
			fid_col  = fid_col + 1
		end

		local pattern = string.rep('[^\t]*\t', fid_col - 1) .. '([^\t]*)'

		for line in lines do
			local facet = string.match(line, pattern)

			local line = remove_nth_element(line, fid_col)

			for i in string.gmatch(facet .. ',', '([^,]*),') do
				if FACETS[i] then
					tables[tonumber(FACETS[i])] = (tables[tonumber(FACETS[i])] or remove_nth_element(header, fid_col)) .. '\\\\' .. line
				end
			end
		end
	else
		for line in lines do
			tables[0] = (tables[0] or header) .. '\\\\' .. line
		end
	end

	return tables
end


-- function to export tables corresponding to the current facet to LaTeX
local function get_facet_table(plot, facet, table_name)
	local data_table = _G['pggf_data_' .. plot]
	local facet_data = data_table[facet] or data_table[0]

	if not facet_data then
		if plot == 'axes' then
			pggf_error('Axis data is misssing for facet ' .. facet .. '. Axis data must be supplied for all facets')
		end
		pggf_warning('No data for plot \'' .. plot .. '\' present in facet ' .. facet, 'emptyfacetwarn')
		tex.sprint(4, '\\pggf@noplotdatatrue')
		return
	end

	if not table_name then
		pggf_error('Missing \'table_name\' \'' .. table_name .. '\' in function \'getTable\'')
	end

	if not string.find(table_name, '^[%a@]+$') then
		pggf_error('Invalid \'table_name\' \'' .. table_name .. '\' in function \'getTable\'')
	end

	tex.sprint('\\pgfplotstableread[format=inline,col sep=tab,row sep=\\\\]{' .. facet_data .. '\\\\}')
	tex.sprint(4, '{\\' .. table_name ..'}')
end


-- function to draw box plots
local function boxplot_from_table(plot_id, facet, plot_type)
	local data_table = _G['pggf_data_boxplot_' .. plot_id]
	local facet_data = data_table[facet] or data_table[0]

	if not facet_data then
		pggf_warning('No data for plot \'' .. plot_id .. '\' present in facet ' .. facet, 'emptyfacetwarn')
		return
	end

	--parse violin data
	local half_violin = token.get_macro('pggf@violin@half')

	local violin_data = {}
	if (plot_type == 'violin') then
		local violin_data_table = _G['pggf_data_violin_' .. plot_id]
		local violin_facet_data = violin_data_table[facet] or violin_data_table[0]

		local lines = string.gmatch(violin_facet_data .. '\\\\', '(.-)\\\\')

		local header = lines()
		local colnames = {}

		for col in string.gmatch(header, '[^\t]*') do
			table.insert(colnames, col)
		end

		for line in lines do
			local line_data = string.gmatch(line, '[^\t]*')
			local line_table = {}

			for i, col in ipairs(colnames) do
				line_table[col] = line_data()
			end

			local sid = line_table['sample_id']
			violin_data[sid] = violin_data[sid] or {}
			violin_data[sid]['x'] = violin_data[sid]['x'] or {}
			violin_data[sid]['y'] = violin_data[sid]['y'] or {}
			table.insert(violin_data[sid]['x'], line_table['x'])
			table.insert(violin_data[sid]['y'], line_table['y'])
		end
	end
	-- end parse violin data

	local shade = token.get_macro('pggf@shade')
	local shade_color
	local shade_percent

	if shade ~= 'false' then
		shade_color = token.get_macro('pggf@shade@color')
		shade_percent = token.get_macro('pggf@shade@percent')
	end

	local lines = string.gmatch(facet_data .. '\\\\', '(.-)\\\\')

	local header = lines()
	local colnames = {}

	for col in string.gmatch(header, '[^\t]*') do
		table.insert(colnames, col)
	end

	local style_col
	local vio_style_col
	if (plot_type == 'violin') then
		vio_style_col = token.get_macro('pggf@stylecol@' .. plot_id)
		style_col = token.get_macro('pggf@boxstylecol@' .. plot_id)
	else
		style_col = token.get_macro('pggf@stylecol@' .. plot_id)
	end

	if style_col then
		style_col = string.gsub(style_col, '%s*([,=])%s*', '%1')
		for col in string.gmatch(style_col, '[^,]*') do
			col = string.gsub(col, '.*=', '')
			if not string.find('\t' .. header .. '\t', '\t' .. col .. '\t') then
				pggf_error('Column \'' .. col .. '\' does not exist.')
			end
		end
	end
	if vio_style_col then
		vio_style_col = string.gsub(vio_style_col, '%s*([,=])%s*', '%1')
		for col in string.gmatch(vio_style_col, '[^,]*') do
			col = string.gsub(col, '.*=', '')
			if not string.find('\t' .. header .. '\t', '\t' .. col .. '\t') then
				pggf_error('Column \'' .. col .. '\' does not exist.')
			end
		end
	end

	local draw_dir
	local orientation
	if string.find('\t' .. header .. '\t', '\tx\t') then
		orientation = 'x'
		draw_dir = 'y'
	else
		draw_dir = 'x'
		orientation = 'y'
	end

	local box_data = {dir = draw_dir}

	local outliers = false
	local points = false

	local samples = token.get_macro('pggf@samples@' .. plot_id)
	local label = token.get_macro('pggf@label@' .. plot_id)
	local mark = token.get_macro('pggf@extramark@' .. plot_id)

	local boxwidth = tonumber(token.get_macro('pggf@boxwidth'))
	local violinwidth = tonumber(token.get_macro('pggf@violinwidth'))
	local halfviolinsep = token.get_macro('pggf@violin@half@sep')

	for line in lines do
		local line_data = string.gmatch(line, '[^\t]*')

		for i, col in ipairs(colnames) do
			box_data[col] = line_data()
		end

		if (not samples) or string.find(',' .. samples .. ',', ',' .. box_data.sample_id .. ',') then
			local box_pos = box_data[orientation]

			if (plot_type == 	'violin') then
				if half_violin == 'true' then
					if tonumber(box_data.box_shift) > 0 then
						box_pos = box_pos + 0.5 * halfviolinsep * violinwidth
					else
						box_pos = box_pos - 0.5 * halfviolinsep * violinwidth
					end
				else
					box_pos = box_pos + box_data.box_shift * violinwidth
				end
			else
				box_pos = box_pos + box_data.box_shift * boxwidth
			end

			-- draw violin plots
			if (plot_type == 'violin') then
				local vio_style = ''
				if vio_style_col then
					for vio_style_col in string.gmatch(vio_style_col, '[^,]*') do
						local opt, col = string.match(vio_style_col, '(.+)=(.+)')
						if opt then
							vio_style = vio_style .. ',' .. opt .. '=' .. box_data[col]
						else
							vio_style = vio_style .. ',' .. box_data[vio_style_col]
						end
					end
					vio_style = string.gsub(vio_style, '^,', '')
				end

				local violin_coords = {}
				local violin_legend = 'violin legend=' .. box_data.dir

				if half_violin == 'true' then
					local halfviolinwidth = (1 - halfviolinsep) * violinwidth
					local lasti
					for i, v in ipairs(violin_data[box_data.sample_id][orientation]) do
						lasti = i
						if i == 1 then
							table.insert(violin_coords, box_pos .. '\t' .. violin_data[box_data.sample_id][draw_dir][i])
						end
						if tonumber(box_data.box_shift) > 0 then
							local coord = box_pos + v * halfviolinwidth .. '\t' .. violin_data[box_data.sample_id][draw_dir][i]
							table.insert(violin_coords, coord)
							violin_legend = 'half violin right legend=' .. box_data.dir
						else
							local coord = box_pos - v * halfviolinwidth .. '\t' .. violin_data[box_data.sample_id][draw_dir][i]
							table.insert(violin_coords, coord)
							violin_legend = 'half violin left legend=' .. box_data.dir
						end
					end
					table.insert(violin_coords, box_pos .. '\t' .. violin_data[box_data.sample_id][draw_dir][lasti])
				else
					for i, v in ipairs(violin_data[box_data.sample_id][orientation]) do
						local coord = box_pos + v * box_data.box_width * violinwidth .. '\t' .. violin_data[box_data.sample_id][draw_dir][i]
						table.insert(violin_coords, i, coord)
						local coord = box_pos - v * box_data.box_width * violinwidth .. '\t' .. violin_data[box_data.sample_id][draw_dir][i]
						table.insert(violin_coords, i, coord)
					end
				end

				local violin_coords_str = table.concat(violin_coords, '\\\\')

				local shaded_color
				if shade == 'left' then
					shaded_color = 'fill=pgffillcolor!' .. 100 - shade_percent * (1 - (box_data[orientation] - 1) / (box_data.n_samples - 1)) .. '!' .. shade_color
				elseif shade == 'right' then
					shaded_color = 'fill=pgffillcolor!' .. 100 - shade_percent * ((box_data[orientation] - 1) / (box_data.n_samples - 1)) .. '!' .. shade_color
				else
					shaded_color = ''
				end

				tex.sprint(
					4,
					'\\addplot+[mark=none,draw=black,fill,',
					violin_legend .. ',pggf@violin@' .. plot_id .. '/.try,',
					vio_style .. ',' .. shaded_color .. ']',
					'table[format=inline,row sep=\\\\,' .. orientation .. '=x,' .. draw_dir .. '=y]{x\ty\\\\',
					violin_coords_str,
					'\\\\} -- cycle;'
				)
			end

			-- draw boxplots
			local hide_boxplots = token.get_macro('pggf@hide@boxplots@' .. plot_id)

			if hide_boxplots and hide_boxplots == 'true' then
				-- do not draw boxplots
			else
				local shaded_color
				local shade_boxplot = 'true'
				if plot_type == 'violin' then
					shade_boxplot = token.get_macro('pggf@shade@boxplot')
				end

				if shade_boxplot == 'true' and shade == 'left' then
					shaded_color = 'fill=pgffillcolor!' .. 100 - shade_percent * (1 - (box_data[orientation] - 1) / (box_data.n_samples - 1)) .. '!' .. shade_color
				elseif shade_boxplot == 'true' and shade == 'right' then
					shaded_color = 'fill=pgffillcolor!' .. 100 - shade_percent * ((box_data[orientation] - 1) / (box_data.n_samples - 1)) .. '!' .. shade_color
				else
					shaded_color = ''
				end

				local style = ''
				if style_col then
					for style_col in string.gmatch(style_col, '[^,]*') do
						local opt, col = string.match(style_col, '(.+)=(.+)')
						if opt then
							style = style .. ',' .. opt .. '=' .. box_data[col]
						else
							style = style .. ',' .. box_data[style_col]
						end
					end
					style = string.gsub(style, '^,', '')
				end

				local forgetplot = ''
				if plot_type == 'violin' then
					forgetplot = 'forget plot'
				end

				local boxplotanchor = 'draw relative anchor=0.5'
				if half_violin == 'true' then
					if tonumber(box_data.box_shift) > 0 then
						boxplotanchor = 'draw relative anchor=0'
					else
						boxplotanchor = 'draw relative anchor=1'
					end
				end

				tex.sprint(
					4,
					'\\addplot+[boxplot prepared={',
					boxplotanchor .. ',',
					'draw direction=' .. box_data.dir .. ',',
					'draw position=' .. box_pos .. ',',
					'box extend=' .. box_data.box_width * boxwidth .. ',',
					'lower whisker=' .. box_data.lw .. ',',
					'lower quartile=' .. box_data.lq .. ',',
					'median=' .. box_data.med .. ',',
					'average=' .. box_data.av .. ',',
					'upper quartile=' .. box_data.uq .. ',',
					'upper whisker=' .. box_data.uw .. '},',
					forgetplot .. ',',
					'boxplot legend,mark=none,draw=black,fill,',
					style .. ',pggf@boxplot@' .. plot_id .. '/.try,' .. shaded_color .. '] coordinates {};'
				)

				-- add coordinate for box maximum
				local coord_pos = box_pos

				if half_violin == 'true' then
					if tonumber(box_data.box_shift) > 0 then
						coord_pos = coord_pos + 0.25 * (halfviolinsep + violinwidth)
					else
						coord_pos = coord_pos - 0.25 * (halfviolinsep + violinwidth)
					end
				end

				if orientation == 'x' then
					coord_pos = coord_pos .. ',' .. box_data.max
				else
					coord_pos = box_data.max .. ',' .. coord_pos
				end

				tex.sprint(4, '\\coordinate (pggf sample ' .. box_data.sample_id .. ') at (' .. coord_pos .. ');')
			end

			-- add sample size
			local sspos = '(' .. box_pos .. ',\\pgfkeysvalueof{/pgfplots/ymin})'
			local ssanc = 'south'
			local ssside = ''

			if orientation == 'y' then
				sspos = '(\\pgfkeysvalueof{/pgfplots/xmin},' .. box_pos .. ')'
				ssanc = 'west'
			end

			if half_violin == 'true' then
				if tonumber(box_data.box_shift) > 0 then
					ssside = ',pggf@samplesize@right/.try'
					ssanc = 'south west'
					if orientation == 'y' then
						ssanc = 'north west'
					end
				else
					ssside = ',pggf@samplesize@left/.try'
					ssanc = 'south east'
					if orientation == 'y' then
						ssanc = 'south west'
					end
				end
			end

			tex.sprint(
				4,
				'\\ifpggf@samplesize',
				'\\node[anchor=' .. ssanc .. ',/pgf/number format/fixed,pggf@samplesize' .. ssside .. '] at ' .. sspos .. ' {\\pgfmathprintnumber{' .. box_data.n .. '}};',
				'\\fi'
			)

			-- add labels
			if label then
				local label_id = 1
				for label in string.gmatch(label, '[^;]*') do
					local label_pos, label_col = string.match(label, '(.+)|(.+)')

					local label_col, symbolic = string.gsub(label_col, '^*', '')

					local label_text = box_data[label_col]
					if not label_text then
						pggf_error('Cannot find column \'' .. label_col .. '\' (in `label(*) = {...}{' .. label_pos .. '}{' .. label_col .. '}`)')
					end

					if symbolic == 0 then
						label_text = '\\pgfmathprintnumber{' .. label_text .. '}'
					end

					if label_pos == 'axis max' then label_pos = '\\pgfkeysvalueof{/pgfplots/' .. draw_dir .. 'max}' end
					if label_pos == 'axis min' then label_pos = '\\pgfkeysvalueof{/pgfplots/' .. draw_dir .. 'min}' end

					if orientation == 'x' then
						label_pos = '(' .. box_pos .. ',' .. (box_data[label_pos] or label_pos) .. ')'
					else
						label_pos = '(' .. (box_data[label_pos] or label_pos) .. ',' .. box_pos .. ')'
					end

					tex.sprint(4, '\\node[/pggf/plot/label_defaults/.try,/pggf/plot/label@' .. plot_id .. '@' .. label_id .. '=' .. box_data.sample_id .. '] at ' .. label_pos .. ' {\\pggf@label@macro{' .. label_text .. '}};')

					label_id = label_id + 1
				end
			end

			-- check if outliers need to be drawn
			if (not outliers and plot_type ~= 'violin') then
				if (tonumber(box_data.outliers) == -1) then
					outliers = true
					points = true
				elseif (tonumber(box_data.outliers) > 0) then
					outliers = true
				end
			end
		end
	end

	-- draw outliers
	if outliers then
		local hide_outliers = token.get_macro('pggf@hide@outliers@' .. plot_id)

		if hide_outliers and hide_outliers == 'true' then
			-- do not add outliers
		else
			local outlier_file = token.get_macro('pggf@data@' .. plot_id) .. '_outlier.tsv'

			if not file_exists(outlier_file) then
				pggf_error('Cannot find outlier file: \'' .. outlier_file .. '\'')
			end

			local outlier_data_tables = split_table_by_facets(outlier_file)
			local outlier_data = outlier_data_tables[facet] or outlier_data_tables[0]

			if samples then
				lines = string.gmatch(outlier_data .. '\\\\', '.-\\\\')
				local data_filtered = lines()
				for line in lines do
					for sample in string.gmatch(samples, '[^,]+') do
						if string.find(line, '\t' .. sample .. '\\\\') then
							data_filtered = data_filtered .. line
						end
					end
				end
				outlier_data = data_filtered
			end

			local jitter = 0
			if points then
				jitter = token.get_macro('pggf@jitter')
			end

			tex.sprint(
				4,
				'\\addplot[only marks,mark=solido,forget plot,pggf@outlier@' .. plot_id .. '/.try]',
				'table[format=inline,row sep=\\\\,'..	orientation .. ' expr=\\thisrow{' .. orientation .. '}+\\thisrow{box_shift}*\\pggf@boxwidth + rand * ' .. jitter .. ',',
				draw_dir .. '=' .. draw_dir .. ']',
				'{' .. outlier_data .. '\\\\};'
			)
		end
	end

	-- add extra marks
	if mark then
		local mark_id = 1
		for mark_col in string.gmatch(mark, '[^;]*') do

			tex.sprint(
				4,
				'\\addplot[only marks,mark=solido,/pggf/plot/extramark@' .. plot_id .. '@' .. mark_id .. ']',
				'table[format=inline,row sep=\\\\,'..	orientation .. ' expr=\\thisrow{' .. orientation .. '}+\\thisrow{box_shift}*\\pggf@boxwidth,',
				draw_dir .. '=' .. mark_col .. ']',
				'{' .. facet_data .. '\\\\};'
			)

			mark_id = mark_id + 1
		end
	end

end


-- function to draw straight lines
local function lines_from_table(plot_id, facet, line_type)
	local data_table = _G['pggf_data_' .. line_type .. '_' .. plot_id]
	local facet_data = data_table[facet] or data_table[0]

	if not facet_data then
		pggf_warning('No data for plot \'' .. plot_id .. '\' present in facet ' .. facet, 'emptyfacetwarn')
		return
	end

	local lines = string.gmatch(facet_data .. '\\\\', '(.-)\\\\')

	local header = lines()
	local colnames = {}

	for col in string.gmatch(header, '[^\t]*') do
		table.insert(colnames, col)
	end

	local style_col

	style_col = token.get_macro('pggf@stylecol@' .. plot_id)

	if style_col then
		style_col = string.gsub(style_col, '%s*([,=])%s*', '%1')
		for col in string.gmatch(style_col, '[^,]*') do
			col = string.gsub(col, '.*=', '')
			if not string.find('\t' .. header .. '\t', '\t' .. col .. '\t') then
				pggf_error('Column \'' .. col .. '\' does not exist.')
			end
		end
	end

	for line in lines do
		local line_data = string.gmatch(line, '[^\t]*')
		local line_table = {}

		for i, col in ipairs(colnames) do
			line_table[col] = line_data()
		end

		if (line_type == 'hline') then
			line_table['slope'] = 0
		end

		local style = ''
		if style_col then
			for style_col in string.gmatch(style_col, '[^,]*') do
				local opt, col = string.match(style_col, '(.+)=(.+)')
				if opt then
					style = style .. ',' .. opt .. '=' .. line_table[col]
				else
					style = style .. ',' .. line_table[style_col]
				end
			end
			style = string.gsub(style, '^,', '')
		end

		if (line_type == 'vline') then
			tex.sprint(
				4,
				'\\draw[pggfplotcolor,/pgfplots/pggf@' .. line_type .. '@' .. plot_id .. '/.try,'.. style .. ']',
				'(' .. line_table.intercept .. ',\\pggfymin) -- (' .. line_table.intercept .. ',\\pggfymax);'
			)
		else
			tex.sprint(
				4,
				'\\draw[pggfplotcolor,/pgfplots/pggf@' .. line_type .. '@' .. plot_id .. '/.try,' .. style .. ']',
				'(\\pggfxmin,' .. line_table.slope .. '*\\pggfxmin+' .. line_table.intercept .. ') -- (\\pggfxmax,' .. line_table.slope .. '*\\pggfxmax+' .. line_table.intercept .. ');'
			)
		end
	end
end

-- function to draw bar plots
local function barplot_from_table(plot_id, facet)
	local data_table = _G['pggf_data_bar_' .. plot_id]
	local facet_data = data_table[facet] or data_table[0]

	if not facet_data then
		pggf_warning('No data for plot \'' .. plot_id .. '\' present in facet ' .. facet, 'emptyfacetwarn')
		return
	end

	local lines = string.gmatch(facet_data .. '\\\\', '(.-)\\\\')

	local header = lines()

	local style_col
	local styles = {}

	style_col = token.get_macro('pggf@stylecol@' .. plot_id)

	if style_col then
		style_col = string.gsub(style_col, '%s*([,=])%s*', '%1')

		for style_col in string.gmatch(style_col, '[^,]*') do
			local opt, col = string.match(style_col, '(.+)=(.+)')
			if not opt then
				col = style_col
			end

			if not string.find('\t' .. header .. '\t', '\t' .. col .. '\t') then
				pggf_error('Column \'' .. col .. '\' does not exist.')
			end

			styles[col] = {}
			styles[col][0] = opt

			local colID = 1
			for colName in string.gmatch(header, '[^\t]*') do
				if colName == col then
					styles[col][-1] = string.rep('[^\t]*\t', colID - 1) .. '([^\t]*)'
					break
				end
				colID  = colID + 1
			end
		end
	end

	local barwidth = token.get_macro('pggf@barwidth')

	local orientation
	local width_col
	local shift_col
	local i = 0

	for colName in string.gmatch(header, '[^\t]*') do
		if colName == 'bar_width' then
			width_col = string.rep('[^\t]*\t', i) .. '([^\t]*)'
		elseif colName == 'x_shift' then
			orientation = 'x'
			shift_col = string.rep('[^\t]*\t', i) .. '([^\t]*)'
		elseif colName == 'y_shift' then
			orientation = 'y'
			shift_col = string.rep('[^\t]*\t', i) .. '([^\t]*)'
		end
		i  = i + 1
	end

	local grouped = false
	local group_col

	if string.find('\t' .. header .. '\t', '\tgroup\t') then
		grouped = true

		local i = 0
		for colName in string.gmatch(header, '[^\t]*') do
			if colName == 'group' then
				group_col = string.rep('[^\t]*\t', i) .. '([^\t]*)'
				break
			end
			i  = i + 1
		end
	end

	local groups = {}
	local group_order = {}
	local bar_data = {}

	for line in lines do
		local group

		if grouped then
			group = string.match(line, group_col)
		else
			group = 'no group'
		end

		if not groups[group] then
			for col in pairs(styles) do
				styles[col][group] = string.match(line, styles[col][-1])
			end

			table.insert(group_order, group)

			groups[group] = {}
			groups[group]['width'] = string.match(line, width_col) * barwidth
			groups[group]['shift'] = string.match(line, shift_col) * barwidth
		end

		bar_data[group] = (bar_data[group] or header) .. '\\\\' .. line
	end

	for _, group in ipairs(group_order) do
		local style = ''
		for col in pairs(styles) do
			if styles[col][0] then
				style = style .. ',' .. styles[col][0] .. '=' .. styles[col][group]
			else
				style = style .. ',' .. styles[col][group]
			end
		end
		style = string.gsub(style, '^,', '')

		local columns

		if orientation == 'x' then
			columns = 'x expr=\\thisrow{x}+' .. groups[group]['shift'] .. ',y=y'
		else
			columns = 'x=x,y expr=\\thisrow{y}+' .. groups[group]['shift']
		end

		tex.sprint(4,
			'\\addplot+[no markers,fill,draw,bar width=' .. groups[group]['width'] .. ',' .. style .. ',pggf@bar@' .. plot_id .. '/.try]',
			'table[format=inline,row sep=\\\\,' .. columns .. ']',
			'{' .. bar_data[group] .. '\\\\};'
		)

		if group ~= 'no group' then
		tex.sprint(4,
			'\\addlegendentry{' .. group .. '}'
		)
		end
	end
end


-- function to draw line plots
local function lineplot_from_table(plot_id, facet)
	local data_table = _G['pggf_data_line_' .. plot_id]
	local facet_data = data_table[facet] or data_table[0]

	if not facet_data then
		pggf_warning('No data for plot \'' .. plot_id .. '\' present in facet ' .. facet, 'emptyfacetwarn')
		return
	end

	local lines = string.gmatch(facet_data .. '\\\\', '(.-)\\\\')

	local header = lines()

	local tables = {}

	local groups = token.get_macro('pggf@line@names')

	local group_col = token.get_macro('pggf@line@grouped@' .. plot_id)

	if group_col then
		if not groups or groups == '' then
			pggf_error('Line plot styles must be created first with `line styles={style name1[style1],style name2, ...}` for grouped line plots to work.')
		end

		if not string.find('\t' .. header .. '\t', '\t' .. group_col .. '\t') then
			pggf_error('Column \'' .. group_col .. '\' does not exist.')
		end

		local group_col_i = 1
		for colName in string.gmatch(header, '[^\t]*') do
			if colName == group_col then
				break
			end
			group_col_i  = group_col_i + 1
		end

		local pattern = string.rep('[^\t]*\t', group_col_i - 1) .. '([^\t]*)'

		for line in lines do
			local group = string.match(line, pattern)

			local line = remove_nth_element(line, group_col_i)

			tables[group] = (tables[group] or remove_nth_element(header, group_col_i)) .. '\\\\' .. line
		end
	else
		tables['no_group'] = facet_data

		if groups and groups ~= '' then
			groups = groups .. ',no_group'
		else
			groups = 'no_group'
		end
	end

	local show_error = token.get_macro('pggf@error@' .. plot_id)
	local error_style = ''
	local error_cols
	local error_dir

	if show_error then
		if not string.find(',lines,ribbon,ribbon and lines,bars,', ',' .. show_error .. ',') then
			pggf_error('Choice \'' .. show_error .. '\' unknown in key \'/pggf/plot/error\'. Use one of \'lines\', \'ribbon\', \'ribbon and lines\', or \'bars\'')
		end

		if show_error == 'bars' then
			error_dir = string.match('\t' .. header .. '\t', '\terror_([xy])_1\t')

			if error_dir then
				error_style = '/pggf/plot/' .. error_dir .. ' error*={error_' .. error_dir .. '_1}{error_' .. error_dir .. '_2}'
			else
				error_dir = string.match('\t' .. header .. '\t', '\terror_([xy])\t')

				if error_dir then
					error_style = '/pggf/plot/' .. error_dir .. ' error={error_' .. error_dir .. '}'
				else
					pggf_error('Can\'t find columns with error data (\'error_' .. error_dir .. '\' or \'error_' .. error_dir .. '_1\' and \'error_' .. error_dir .. '_2\')')
				end
			end
		else
			error_dir = string.match('\t' .. header .. '\t', '\terror_([xy])_min\t')

			if not error_dir or not string.find('\t' .. header .. '\t', '\terror_' .. error_dir .. '_max\t') then
				pggf_error('Can\'t find columns with error data (\'error_' .. error_dir .. '_min\' and \'error_' .. error_dir .. '_max\')')
			end

			if error_dir == 'y' then
				error_cols = 'x=x,y=error_y_'
			else
				error_cols = 'y=y,x=error_x_'
			end
		end
	end

	groups = string.gmatch(groups .. ',', '(.-),')

	for group in groups do
		local group_style

		if group == 'no_group' then
			group_style = ''
		else
			group_style = group
		end

		if tables[group] then
			if show_error == 'ribbon' or show_error == 'ribbon and lines' then
				local lines = string.gmatch(tables[group] .. '\\\\', '(.-)\\\\')

				local header = lines()

				local i = 1
				local min_col_i
				local max_col_i

				for colName in string.gmatch(header, '[^\t]*') do
					if colName == 'error_' .. error_dir .. '_min' then
						min_col_i = i
					elseif colName == 'error_' .. error_dir .. '_max' then
						max_col_i = i
					end
					i  = i + 1
				end

				local error_table = {}

				local i = 1

				for line in lines do
					table.insert(error_table, i, remove_nth_element(line, max_col_i))
					line = switch_elements(line, min_col_i, max_col_i)
					table.insert(error_table, i, remove_nth_element(line, max_col_i))
					i = i + 1
				end

				local error_data

				for i, line in ipairs(error_table) do
					error_data = (error_data or remove_nth_element(header, max_col_i)) .. '\\\\' .. line
				end

				tex.sprint(
					4,
					'\\addplot[mark=none,draw=none,fill=pggfplotcolor,' .. group_style .. ',pggf@line@' .. plot_id .. '/.try,',
					'forget plot,fill opacity=0.2,pggf@line@error@ribbon@' .. plot_id ..'/.try]',
					'table[format=inline,row sep=\\\\,' .. error_cols .. 'min] {' .. error_data ..'\\\\} -- cycle;'
				)
			end

			if show_error == 'lines' or show_error == 'ribbon and lines' then
				tex.sprint(
					4,
					'\\addplot[mark=none,pggfplotcolor,' .. group_style .. ',pggf@line@' .. plot_id .. '/.try,',
					'forget plot,dashed,pggf@line@error@min@' .. plot_id ..'/.try]',
					'table[format=inline,row sep=\\\\,' .. error_cols .. 'min] {' .. tables[group] ..'\\\\};'
				)
				tex.sprint(
					4,
					'\\addplot[mark=none,pggfplotcolor,' .. group_style .. ',pggf@line@' .. plot_id .. '/.try,',
					'forget plot,dashed,pggf@line@error@max@' .. plot_id ..'/.try]',
					'table[format=inline,row sep=\\\\,' .. error_cols .. 'max] {' .. tables[group] ..'\\\\};'
				)
			end

			tex.sprint(
				4,
				'\\addplot[mark=none,pggfplotcolor,' .. group_style .. ',' .. error_style .. ',pggf@line@' .. plot_id .. '/.try]',
				'table[format=inline,row sep=\\\\,x=x,y=y] {' .. tables[group] ..'\\\\};'
			)
			tables[group] = nil
		else
			tex.sprint('\\addlegendimage{' .. group .. '}')
		end
	end

	for group in pairs(tables) do
		pggf_warning('Group \'' .. group .. '\' is found in the data but not in `line styles`. It will not be plotted. Use `missing group warning = false` to silence this warning.', 'missinggroupwarn')
	end
end


-- function to draw density plots
local function density_from_table(plot_id, facet)
	local data_table = _G['pggf_data_density_' .. plot_id]
	local facet_data = data_table[facet] or data_table[0]

	if not facet_data then
		pggf_warning('No data for plot \'' .. plot_id .. '\' present in facet ' .. facet, 'emptyfacetwarn')
		return
	end

	local lines = string.gmatch(facet_data .. '\\\\', '(.-)\\\\')

	local header = lines()

	local tables = {}

	local groups = token.get_macro('pggf@density@names')

	local group_col = token.get_macro('pggf@density@grouped@' .. plot_id)

	if group_col then
		if not groups or groups == '' then
			pggf_error('Density plot styles must be created first with `density styles={style name1[style1],style name2, ...}` for grouped density plots to work.')
		end

		if not string.find('\t' .. header .. '\t', '\t' .. group_col .. '\t') then
			pggf_error('Column \'' .. group_col .. '\' does not exist.')
		end

		local group_col_i = 1
		for colName in string.gmatch(header, '[^\t]*') do
			if colName == group_col then
				break
			end
			group_col_i  = group_col_i + 1
		end

		local pattern = string.rep('[^\t]*\t', group_col_i - 1) .. '([^\t]*)'

		for line in lines do
			local group = string.match(line, pattern)

			local line = remove_nth_element(line, group_col_i)

			tables[group] = (tables[group] or remove_nth_element(header, group_col_i)) .. '\\\\' .. line
		end
	else
		tables['no_group'] = facet_data

		if groups and groups ~= '' then
			groups = groups .. ',no_group'
		else
			groups = 'no_group'
		end
	end

	groups = string.gmatch(groups .. ',', '(.-),')

	for group in groups do
		local group_style

		if group == 'no_group' then
			group_style = ''
		else
			group_style = group
		end

		if tables[group] then
			tex.sprint(
				4,
				'\\addplot+[mark=none,draw,fill,' .. group_style .. ',pggf@density@' .. plot_id .. '/.try]',
				'table[format=inline,row sep=\\\\,x=x,y=y] {' .. tables[group] ..'\\\\}\\ifpggf@close@density\\closedcycle\\fi;'
			)
			tables[group] = nil
		else
			tex.sprint('\\addlegendimage{' .. group .. '}')
		end
	end

	for group in pairs(tables) do
		pggf_warning('Group \'' .. group .. '\' is found in the data but not in `density styles`. It will not be plotted. Use `missing group warning = false` to silence this warning.', 'missinggroupwarn')
	end
end


-- function to draw ridgeline plots
local function ridgeline_from_table(plot_id, facet)
	local data_table = _G['pggf_data_ridgeline_' .. plot_id]
	local facet_data = data_table[facet] or data_table[0]

	if not facet_data then
		pggf_warning('No data for plot \'' .. plot_id .. '\' present in facet ' .. facet, 'emptyfacetwarn')
		return
	end

	local lines = string.gmatch(facet_data .. '\\\\', '(.-)\\\\')

	local header = lines()

	local style_col
	local styles = {}

	style_col = token.get_macro('pggf@stylecol@' .. plot_id)

	if style_col then
		style_col = string.gsub(style_col, '%s*([,=])%s*', '%1')

		for style_col in string.gmatch(style_col, '[^,]*') do
			local opt, col = string.match(style_col, '(.+)=(.+)')
			if not opt then
				col = style_col
			end

			if not string.find('\t' .. header .. '\t', '\t' .. col .. '\t') then
				pggf_error('Column \'' .. col .. '\' does not exist.')
			end

			styles[col] = {}
			styles[col][0] = opt

			local colID = 1
			for colName in string.gmatch(header, '[^\t]*') do
				if colName == col then
					styles[col][-1] = string.rep('[^\t]*\t', colID - 1) .. '([^\t]*)'
					break
				end
				colID  = colID + 1
			end
		end
	end

	local tables = {}

	local group_col = ''
	local mapping = ''

	if string.find('\t' .. header .. '\t', '\tdensity_y\t') then
		group_col = 'y'
		mapping = 'x=x,y=density_y'
	elseif string.find('\t' .. header .. '\t', '\tdensity_x\t') then
		group_col = 'x'
		mapping = 'x=density_x,y=y'
	else
		pggf_error('Column \'density_y\' or \'density_x\' is required but missing.')
	end

	local group_col_i = 1
	for colName in string.gmatch(header, '[^\t]*') do
		if colName == group_col then
			break
		end
		group_col_i  = group_col_i + 1
	end

	local group_pattern = string.rep('[^\t]*\t', group_col_i - 1) .. '([^\t]*)'

	local density_col_i = 1
	for colName in string.gmatch(header, '[^\t]*') do
		if colName == 'density_' .. group_col then
			break
		end
		density_col_i  = density_col_i + 1
	end

	local density_pattern = string.rep('[^\t]*\t', density_col_i - 1) .. '([^\t]*)'

	local scale_factor = token.get_macro('pggfridgescale')

	local group_order = {}

	for line in lines do
		local group = string.match(line, group_pattern)

		if not tables[group] then
			table.insert(group_order, group)

			for col in pairs(styles) do
				styles[col][group] = string.match(line, styles[col][-1])
			end
		end

		local line = modify_nth_element(line, density_col_i, scale_factor, group)

		tables[group] = (tables[group] or header) .. '\\\\' .. line
	end

	for _, group in ipairs(group_order) do
		local style = ''
		for col in pairs(styles) do
			if styles[col][0] then
				style = style .. ',' .. styles[col][0] .. '=' .. styles[col][group]
			else
				style = style .. ',' .. styles[col][group]
			end
		end
		style = string.gsub(style, '^,', '')

		tex.sprint(
			4,
			'\\addplot+[mark=none,draw,fill,' .. style .. ',pggf@ridgeline@' .. plot_id .. '/.try]',
			'table[format=inline,row sep=\\\\,' .. mapping .. '] {' .. tables[group] ..'\\\\}\\ifpggf@close@density--cycle\\fi;'
		)
	end
end


-- helper function to find the next free level for a significance bracket
local function find_level(s1, s2, used_levels)
	local searching = true
	local level = 0

	if token.get_macro('pggf@signif@shiftsame') == 'false' then
		if s2 - s1 == 1 then
			level = 1
			searching = false
		else
			s1 = s1 + 1
			s2 = s2 -1
		end
	end

	while searching do
		local i = s1
		level = level + 1
		while i <= s2 and (not used_levels[i] or not used_levels[i][level]) do
			if i == s2 then
				searching = false
			end
			i = i + 1
		end
	end

	return(level)
end


-- function to draw significance brackets
local function signif_from_table(plot_id, facet)
	local data_table = _G['pggf_data_signif_' .. plot_id]
	local facet_data = data_table[facet] or data_table[0]

	if not facet_data then
		pggf_warning('No data for significance brackets \'' .. plot_id .. '\' present in facet ' .. facet, 'emptyfacetwarn')
		return
	end

	local lines = string.gmatch(facet_data .. '\\\\', '(.-)\\\\')

	local header = lines()
	local colnames = {}

	for col in string.gmatch(header, '[^\t]*') do
		table.insert(colnames, col)
	end

	local orientation
	local draw_dir
	if string.find('\t' .. header .. '\t', '\tmax_y\t') then
		orientation = 'x'
		draw_dir = 'y'
	else
		orientation = 'y'
		draw_dir = 'x'
	end

	local signif_data = {}

	local hide_brackets = token.get_macro('pggf@hide@brackets@' .. plot_id)

	local manual_levels = token.get_macro('pggf@signif@levels@' .. plot_id)

	local level_table = {}

	for line in lines do
		local line_data = string.gmatch(line, '[^\t]*')

		for i, col in ipairs(colnames) do
			signif_data[col] = line_data()
		end

		local cur_bracket = signif_data.sample1 .. '-' .. signif_data.sample2

		local show_bracket = true
		if hide_brackets and string.find(',' .. hide_brackets .. ',', ',' .. cur_bracket .. ',', 1, true) then
				show_bracket = false
		end

		-- test if level was set manually
		local level
		local manual_level = false

		if manual_levels and string.find(',' .. manual_levels, ',' .. cur_bracket .. '=', 1, true) then
			for level_spec in string.gmatch(manual_levels, '[^\\,]*') do
				if string.match(level_spec, cur_bracket .. '=(.*)') then
					level = string.match(level_spec, cur_bracket .. '=(.*)')
					manual_level = true
				end
			end
		end

		-- find next free level for significance bracket (if not ignoring hidden brackets)
		if not manual_level then
			if token.get_macro('pggf@signif@ignorehidden') == 'false' then
				level = find_level(tonumber(signif_data.sample1), tonumber(signif_data.sample2), level_table)

				for i = tonumber(signif_data.sample1), tonumber(signif_data.sample2) do
					level_table[i] = level_table[i] or {}
					level_table[i][level] = 1
				end
			end
		end

		if show_bracket then
			-- find next free level for significance bracket (ignoring hidden brackets)
			if not manual_level then
				if token.get_macro('pggf@signif@ignorehidden') == 'true' then
					level = find_level(tonumber(signif_data.sample1), tonumber(signif_data.sample2), level_table)

					for i = tonumber(signif_data.sample1), tonumber(signif_data.sample2) do
						level_table[i] = level_table[i] or {}
						level_table[i][level] = 1
					end
				end
			end

			-- determine level coordinates
			local level_coords

			if orientation == 'x' then
				level_coords = 'axis direction cs:0,{' .. level - 1 .. '\\pggf@signif@leveldist/\\pgfplotsunitylength/transformdirectiony(1)}'
			else
				level_coords = 'axis direction cs:{' .. level - 1 .. '\\pggf@signif@leveldist/\\pgfplotsunitxlength/transformdirectionx(1)}, 0'
			end

			-- determine how much to raise the significance brackets
			local raise = string.gsub(token.get_macro('pggf@signif@raise'), '%s*=%s*', '=')
			local raise_type, raise_value = string.match(raise, '(.+)=(.+)')
			local raise_coords = 'axis direction cs:'

			if raise_type == 'rel' then
				if tonumber(raise_value) then
					raise_value = tonumber(raise_value) * tonumber(signif_data.range)
				else
					pggf_error('The argument of `raise = {' .. raise_type .. ' = ...}` must be a number')
				end
			elseif raise_type =='abs' then
				if not tonumber(raise_value) then
					pggf_error('The argument of `raise = {' .. raise_type .. ' = ...}` must be a number')
				end
			elseif raise_type == 'len' then
				raise_value = '{' .. raise_value .. '/\\pgfplotsunit' .. draw_dir .. 'length/transformdirection' .. draw_dir .. '(1)}'
			else
				pggf_error('The argument of `raise` must be one of: "rel = <number>", "abs = <number>", or "len = <length>"')
			end

			if orientation == 'x' then
				raise_coords = raise_coords .. '0,' .. raise_value
			else
				raise_coords = raise_coords .. raise_value .. ',0'
			end

			-- determine how much the significance brackets extend below the top line
			local extend = string.gsub(token.get_macro('pggf@signif@extend'), '%s*=%s*', '=')
			local extend_coords

			if string.find(extend, '^%s*sample%s*$') then
				extend_coords = ''
			else
				extend_coords = '|- 1,' .. signif_data['max_' .. draw_dir] .. ') + (' .. level_coords .. ') + ('
				local extend_type, extend_value = string.match(extend, '(.+)=(.+)')
				if extend_type == 'rel' then
					if tonumber(extend_value) then
						extend_value = tonumber(extend_value) * tonumber(signif_data.range)
					else
						pggf_error('The argument of `extend = {' .. extend_type .. ' = ...}` must be a number')
					end
				elseif extend_type =='abs' then
					if not tonumber(extend_value) then
						pggf_error('The argument of `extend = {' .. extend_type .. ' = ...}` must be a number')
					end
				elseif extend_type == 'len' then
					extend_value = '{' .. extend_value .. '/\\pgfplotsunit' .. draw_dir .. 'length/transformdirection' .. draw_dir .. '(1)}'
				else
					pggf_error('The argument of `extend` must be one of: "sample", "rel = <number>", "abs = <number>", or "len = <length>"')
				end

				if orientation == 'x' then
					extend_coords = extend_coords .. 'axis direction cs:0,-' .. extend_value
				else
					extend_coords = extend_coords .. 'axis direction cs:-' .. extend_value .. ',0'
				end
			end

			-- draw bracket
			tex.sprint(
				4,
				'\\draw[/pgfplots/pggf@signif@' .. plot_id .. '/.try,',
				'/pgfplots/pggf@signif@' .. plot_id .. '@' .. signif_data.sample1 .. '-' .. signif_data.sample2 .. '/.try]',
				'($(pggf sample ' .. signif_data.sample1 .. extend_coords .. ') + (' .. raise_coords ..')$) |-',
				'($(pggf sample ' .. signif_data.sample2 .. ' |- 1,' .. signif_data['max_' .. draw_dir] .. ') + (' .. raise_coords ..') + (' .. level_coords .. ')$)', 
				'node[pos=0.75,anchor=south,/tikz/pggf@signif@label@' .. plot_id .. '/.try] {\\pggf@print@signif{' .. signif_data.p_value .. '}} --',
				'($(pggf sample ' .. signif_data.sample2 .. extend_coords .. ') + (' .. raise_coords ..')$);'
			)
		end
	end
end


-- make functions available in LaTeX
return {
	splitTable = split_table_by_facets,
	getTable = get_facet_table,
	getFacets = get_facet_info,
	boxPlot = boxplot_from_table,
	simpleLines = lines_from_table,
	barPlot = barplot_from_table,
	linePlot = lineplot_from_table,
	densityPlot = density_from_table,
	ridgelinePlot = ridgeline_from_table,
	signifBracket = signif_from_table
}