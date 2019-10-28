'''
Created on 29-July-2019

@author: Darshan
'''
from numpy import math
from Connections.connection_calculations import ConnectionCalculations
import svgwrite
import cairosvg
import numpy as np
import os


class Tension_drawing(object):
	def __init__(self, input_dict, output_dict, memb_data, folder):
		"""

		Args:
			input_dict: input parameters from GUI
			output_dict:  output parameters based on calculation
			beam_data:  geometric properties of beam
			folder: path to save the generated images

		Returns:
			None

		"""
		print "calculation", input_dict
		self.folder = folder
		self.section_type = input_dict["Member"]["SectionType"]
		self.conn_loc = input_dict["Member"]["Location"]
		self.act_member_length = float(input_dict["Member"]["Member_length"])
		self.member_length = 750
		self.weld_inline = float(input_dict["Weld"]["inline_tension"])
		self.weld_oppline = float(input_dict["Weld"]["oppline_tension"])
		self.beam_designation = memb_data['Designation']

		if input_dict["Member"]["SectionType"] == "Angles":
			self.member_leg = memb_data["AXB"]
			self.leg = self.member_leg.split("x")
			self.leg1 = self.leg[0]
			self.leg2 = self.leg[1]
			self.leg_min = min(float(self.leg1),float(self.leg2))
			self.leg_max = max(float(self.leg1),float(self.leg2))
			self.t = float(memb_data["t"])

		else:
			self.member_tw = float(memb_data["tw"])
			self.member_tf = float(memb_data["T"])
			self.member_d = float(memb_data["D"])
			self.member_B = float(memb_data["B"])

		if self.conn_loc == "Back to Back Angles" or self.conn_loc == "Back to Back Web" or self.conn_loc == "Star Angles":
			self.plate_thickness = float(input_dict["Weld"]["Platethickness"])


	def add_s_marker(self, dwg):
		"""

		Args:
			dwg: svgwrite (obj)

		Returns:
			Container for all svg elements

		"""
		smarker = dwg.marker(insert=(8, 3), size=(30, 30), orient="auto")
		smarker.add(dwg.path(d=" M0,0 L3,3 L0,8 L8,3 L0,0", fill="black"))
		dwg.defs.add(smarker)
		return smarker

	def add_section_marker(self, dwg):
		"""

		Args:
			dwg: svgwrite (obj)

		Returns:
			Container for all svg elements

		"""
		section_marker = dwg.marker(insert=(0, 5), size=(10, 10), orient="auto")
		section_marker.add(dwg.path(d="M 0 0 L 10 5 L 0 10 z", fill="blue", stroke="black"))
		dwg.defs.add(section_marker)
		return section_marker

	def add_e_marker(self, dwg):
		"""

		Args:
			dwg: svgwrite (obj)

		Returns:
			Container for all svg elements

		"""
		emarker = dwg.marker(insert=(0, 3), size=(30, 20), orient="auto")
		emarker.add(dwg.path(d=" M0,3 L8,8 L5,3 L8,0 L0,3", fill="black"))
		dwg.defs.add(emarker)
		return emarker

	def draw_start_arrow(self, line, s_arrow):
		"""

		Args:
			line: start line marker
			s_arrow: start arrow

		Returns:
			None

		"""
		line["marker-start"] = s_arrow.get_funciri()

	def draw_end_arrow(self, line, e_arrow):
		"""

		Args:
			line: end line marker
			e_arrow: end arrow

		Returns:
			None

		"""
		line["marker-end"] = e_arrow.get_funciri()

	def draw_faint_line(self, pt_one, pt_two, dwg):
		"""

		Args:
			pt_one: first point
			pt_two: second point
			dwg: svgwrite (obj)

		Returns:
			None

		"""
		dwg.add(dwg.line(pt_one, pt_two).stroke("#D8D8D8", width=2.5, linecap="square", opacity=0.70))

	def draw_dimension_outer_arrow(self, dwg, pt1, pt2, text, params):
		# TODO

		"""

		Args:
			dwg: svgwrite (obj)
			pt1: first point
			pt2: second point
			text: text message
			params:

		Returns:
			None

		"""
		smarker = self.add_s_marker(dwg)
		emarker = self.add_e_marker(dwg)
		line_vector = pt2 - pt1  # [a, b]
		normal_vector = np.array([-line_vector[1], line_vector[0]])  # [-b, a]
		normal_unit_vector = self.normalize(normal_vector)

		if params["lineori"] == "left":
			normal_unit_vector = -normal_unit_vector

		Q1 = pt1 + params["offset"] * normal_unit_vector
		Q2 = pt2 + params["offset"] * normal_unit_vector
		line = dwg.add(dwg.line(Q1, Q2).stroke("black", width=2.5, linecap="square"))
		self.draw_start_arrow(line, emarker)
		self.draw_end_arrow(line, smarker)

		Q12_mid = 0.5 * (Q1 + Q2)
		text_pt = Q12_mid + params["textoffset"] * normal_unit_vector - (len(text)*16)/2
		dwg.add(dwg.text(text, insert=text_pt, fill="black", font_family="sans-serif", font_size=28))

		L1 = Q1 + params["endlinedim"] * normal_unit_vector
		L2 = Q1 + params["endlinedim"] * (-normal_unit_vector)
		dwg.add(dwg.line(L1, L2).stroke("black", width=2.5, linecap="square", opacity=1.0))

		L3 = Q2 + params["endlinedim"] * normal_unit_vector
		L4 = Q2 + params["endlinedim"] * (-normal_unit_vector)
		dwg.add(dwg.line(L3, L4).stroke("black", width=2.5, linecap="square", opacity=1.0))

	def draw_dimension_outer_arrow_side(self, dwg, pt1, pt2, text, params):
		# TODO

		"""

		Args:
			dwg: svgwrite (obj)
			pt1: first point
			pt2: second point
			text: text message
			params:

		Returns:
			None

		"""
		smarker = self.add_s_marker(dwg)
		emarker = self.add_e_marker(dwg)
		line_vector = pt2 - pt1  # [a, b]
		normal_vector = np.array([-line_vector[1], line_vector[0]])  # [-b, a]
		normal_unit_vector = self.normalize(normal_vector)

		if params["lineori"] == "left":
			normal_unit_vector = -normal_unit_vector

		Q1 = pt1 + params["offset"] * normal_unit_vector
		Q2 = pt2 + params["offset"] * normal_unit_vector
		line = dwg.add(dwg.line(Q1, Q2).stroke("black", width=1, linecap="square"))
		self.draw_start_arrow(line, emarker)
		self.draw_end_arrow(line, smarker)

		Q12_mid = 0.5 * (Q1 + Q2)
		text_pt = Q12_mid + params["textoffset"] * normal_unit_vector

		if self.section_type == "Angles":
			dwg.add(dwg.text(text, insert=text_pt, fill="black", font_family="sans-serif", font_size=15))
		else:
			dwg.add(dwg.text(text, insert=text_pt, fill="black", font_family="sans-serif", font_size=28))


		L1 = Q1 + params["endlinedim"] * normal_unit_vector
		L2 = Q1 + params["endlinedim"] * (-normal_unit_vector)
		dwg.add(dwg.line(L1, L2).stroke("black", width=1, linecap="square", opacity=1.0))

		L3 = Q2 + params["endlinedim"] * normal_unit_vector
		L4 = Q2 + params["endlinedim"] * (-normal_unit_vector)
		dwg.add(dwg.line(L3, L4).stroke("black", width=1, linecap="square", opacity=1.0))

	def normalize(self, vector):
		"""

		Args:
			vector: list containing X, Y ordinates of vector

		Returns:
			vector containing normalized X and Y ordinates

		"""
		a = vector[0]
		b = vector[1]
		magnitude = math.sqrt(a * a + b * b)
		return vector / magnitude

	def draw_cross_section(self, dwg, pt_a, pt_b, text_pt, text):
		"""

		Args:
			dwg: svgwrite (obj)
			pt_a: point A
			pt_b: point B
			text_pt: text point
			text: text message

		Returns:
			None

		"""
		line = dwg.add(dwg.line(pt_a, pt_b).stroke("black", width=2.5, linecap="square"))
		sec_arrow = self.add_section_marker(dwg)
		self.draw_end_arrow(line, sec_arrow)
		dwg.add(dwg.text(text, insert=text_pt, fill="black", font_family="sans-serif", font_size=40))

	def draw_dimension_inner_arrow(self, dwg, pt_a, pt_b, text, params):
		# TODO
		"""

		Args:
			dwg: svgwrite (obj)
			pt_a: point A
			pt_b: point B
			text: text message
			params:
				params["offset"] (float): offset of the dimension line
				params["textoffset"] (float): offset of text from dimension line
				params["lineori"] (float): orientation of line [right/left]
				params["endlinedim"] (float): dimension line at the end of the outer arrow
				params["arrowlen"] (float): size of the arrow

		Returns:
			None

		"""
		smarker = self.add_s_marker(dwg)
		emarker = self.add_e_marker(dwg)
		u = pt_b - pt_a  # [a, b]
		u_unit = self.normalize(u)
		v_unit = np.array([-u_unit[1], u_unit[0]])  # [-b, a]

		A1 = pt_a + params["endlinedim"] * v_unit
		A2 = pt_a + params["endlinedim"] * (-v_unit)
		dwg.add(dwg.line(A1, A2).stroke("black", width=2.5, linecap="square"))

		B1 = pt_b + params["endlinedim"] * v_unit
		B2 = pt_a + params["endlinedim"] * (-v_unit)
		dwg.add(dwg.line(B1, B2).stroke("black", width=2.5, linecap="square"))

		A3 = pt_a - params["arrowlen"] * u_unit
		B3 = pt_b + params["arrowlen"] * u_unit
		line = dwg.add(dwg.line(A3, pt_a).stroke("black", width=2.5, linecap="square"))
		self.draw_end_arrow(line, smarker)

		line = dwg.add(dwg.line(B3, pt_b).stroke("black", width=2.5, linecap="butt"))
		self.draw_end_arrow(line, smarker)

		if params["lineori"] == "right":
			text_pt = B3 + params["textoffset"] * u_unit
		else:
			text_pt = A3 - (params["textoffset"] + 100) * u_unit

		dwg.add(dwg.text(text, insert=text_pt, fill="black", font_family='sans-serif', font_size=28))

	def draw_oriented_arrow(self, dwg, point, theta, orientation, offset, textup, textdown, element):
		"""

		Args:
			dwg: svgwrite (obj)
			point: point
			theta: theta
			orientation: direction (east, west, south, north)
			offset: position of the text
			textup: text written above line
			textdown: text written below line

		Returns:
			None

		"""
		# Right Up.
		theta = math.radians(theta)
		char_width = 14
		x_vector = np.array([1, 0])
		y_vector = np.array([0, 1])

		p1 = point
		length_A = offset / (math.sin(theta))

		arrow_vector = None
		if orientation == "NE":
			arrow_vector = np.array([-math.cos(theta), math.sin(theta)])
		elif orientation == "NW":
			arrow_vector = np.array([math.cos(theta), math.sin(theta)])
		elif orientation == "SE":
			arrow_vector = np.array([-math.cos(theta), -math.sin(theta)])
		elif orientation == "SW":
			arrow_vector = np.array([math.cos(theta), -math.sin(theta)])
		p2 = p1 - length_A * arrow_vector

		text = textdown if len(textdown) > len(textup) else textup
		length_B = len(text) * char_width

		label_vector = None
		if orientation == "NE":
			label_vector = -x_vector
		elif orientation == "NW":
			label_vector = x_vector
		elif orientation == "SE":
			label_vector = -x_vector
		elif orientation == "SW":
			label_vector = x_vector
		p3 = p2 + length_B * (-label_vector)

		text_offset = 18
		offset_vector = -y_vector

		text_point_up = None
		text_point_down = None
		if orientation == "NE":
			text_point_up = p2 + 0.2 * length_B * (-label_vector) + text_offset * offset_vector
			text_point_down = p2 - 0.2 * length_B * label_vector - (text_offset + 15) * offset_vector
		elif orientation == "NW":
			text_point_up = p3 + 0.05 * length_B * (label_vector) + text_offset * offset_vector
			text_point_down = p3 - 0.05 * length_B * label_vector - (text_offset + 15) * offset_vector
		elif orientation == "SE":
			text_point_up = p2 + 0.2 * length_B * (-label_vector) + text_offset * offset_vector
			text_point_down = p2 - 0.2 * length_B * label_vector - (text_offset + 15) * offset_vector
		elif orientation == "SW":
			text_point_up = p3 + 0.05 * length_B * (label_vector) + text_offset * offset_vector
			text_point_down = p3 - 0.05 * length_B * label_vector - (text_offset + 15) * offset_vector

		line = dwg.add(dwg.polyline(points=[p1, p2, p3], fill="none", stroke='black', stroke_width=2.5))

		emarker = self.add_e_marker(dwg)
		self.draw_start_arrow(line, emarker)

		dwg.add(dwg.text(textup, insert=text_point_up, fill='black', font_family='sans-serif', font_size=28))
		dwg.add(dwg.text(textdown, insert=text_point_down, fill='black', font_family='sans-serif', font_size=28))

		if element == "weld":
			if self.weld == "Fillet Weld":
				if orientation == "NE":
					self.draw_weld_marker1(dwg, 30, 7.5, line)
				else:
					self.draw_weld_marker2(dwg, 30, 7.5, line)
			else:
				if orientation == "NE":
					self.draw_weld_marker3(dwg, 15, -8.5, line)
				else:
					self.draw_weld_marker4(dwg, 15, 8.5, line)

			if self.stiffener_weld == 1:
				if orientation == "NE":
					self.draw_weld_marker1(dwg, 30, 7.5, line)
				else:
					self.draw_weld_marker2(dwg, 30, 7.5, line)
			else:
				pass

			print "successful"

	def draw_weld_marker1(self, dwg, oriX, oriY, line):
		weldMarker = dwg.marker(insert=(oriX, oriY), size=(15, 15), orient="auto")
		# weldMarker.add(dwg.path(d="M 0 0 L 8 7.5 L 0 15 z", fill='none', stroke='black'))
		weldMarker.add(dwg.path(d="M 15 7.5 L 8 0 L 8 15 z", fill='none', stroke='black'))
		dwg.defs.add(weldMarker)
		self.draw_end_arrow(line, weldMarker)

	def draw_weld_marker2(self, dwg, oriX, oriY, line):
		weldMarker = dwg.marker(insert=(oriX, oriY), size=(15, 15), orient="auto")
		# weldMarker.add(dwg.path(d="M 0 0 L 8 7.5 L 0 15 z", fill='none', stroke='black'))
		weldMarker.add(dwg.path(d="M 0 7.5 L 8 0 L 8 15 z", fill='none', stroke='black'))
		dwg.defs.add(weldMarker)
		self.draw_end_arrow(line, weldMarker)

	def draw_weld_marker3(self, dwg, oriX, oriY, line):
		weldMarker = dwg.marker(insert=(oriX, oriY), size=(15, 15), orient="auto")
		# weldMarker.add(dwg.path(d="M 0 0 L 8 7.5 L 0 15 z", fill='none', stroke='black'))
		weldMarker.add(dwg.path(d="M 0 0 L 0 -7.5 L 7.5 0 ", fill='none', stroke='black'))
		dwg.defs.add(weldMarker)
		self.draw_end_arrow(line, weldMarker)

	def draw_weld_marker4(self, dwg, oriX, oriY, line):
		weldMarker = dwg.marker(insert=(oriX, oriY), size=(15, 15), orient="auto")
		# weldMarker.add(dwg.path(d="M 0 0 L 8 7.5 L 0 15 z", fill='none', stroke='black'))
		weldMarker.add(dwg.path(d="M 0 0 L 0 7.5 L -7.5 0 ", fill='none', stroke='black'))
		dwg.defs.add(weldMarker)
		self.draw_end_arrow(line, weldMarker)

	def draw_oriented_arrow_side(self, dwg, point, theta, orientation, offset, textup, textdown, element):
		"""

		Args:
			dwg: svgwrite (obj)
			point: point
			theta: theta
			orientation: direction (east, west, south, north)
			offset: position of the text
			textup: text written above line
			textdown: text written below line

		Returns:
			None

		"""
		# Right Up.
		theta = math.radians(theta)
		char_width = 10
		x_vector = np.array([1, 0])
		y_vector = np.array([0, 1])

		p1 = point
		length_A = offset / (math.sin(theta))

		arrow_vector = None
		if orientation == "NE":
			arrow_vector = np.array([-math.cos(theta), math.sin(theta)])
		elif orientation == "NW":
			arrow_vector = np.array([math.cos(theta), math.sin(theta)])
		elif orientation == "SE":
			arrow_vector = np.array([-math.cos(theta), -math.sin(theta)])
		elif orientation == "SW":
			arrow_vector = np.array([math.cos(theta), -math.sin(theta)])
		p2 = p1 - length_A * arrow_vector

		text = textdown if len(textdown) > len(textup) else textup
		length_B = len(text) * char_width

		label_vector = None
		if orientation == "NE":
			label_vector = -x_vector
		elif orientation == "NW":
			label_vector = x_vector
		elif orientation == "SE":
			label_vector = -x_vector
		elif orientation == "SW":
			label_vector = x_vector
		p3 = p2 + length_B * (-label_vector)

		text_offset = 18
		offset_vector = -y_vector

		text_point_up = None
		text_point_down = None
		if orientation == "NE":
			text_point_up = p2 + 0.2 * length_B * (-label_vector) + text_offset * offset_vector
			text_point_down = p2 - 0.2 * length_B * label_vector - (text_offset + 15) * offset_vector
		elif orientation == "NW":
			text_point_up = p3 + 0.05 * length_B * (label_vector) + text_offset * offset_vector
			text_point_down = p3 - 0.05 * length_B * label_vector - (text_offset + 15) * offset_vector
		elif orientation == "SE":
			text_point_up = p2 + 0.2 * length_B * (-label_vector) + text_offset * offset_vector
			text_point_down = p2 - 0.2 * length_B * label_vector - (text_offset + 15) * offset_vector
		elif orientation == "SW":
			text_point_up = p3 + 0.05 * length_B * (label_vector) + text_offset * offset_vector
			text_point_down = p3 - 0.05 * length_B * label_vector - (text_offset + 15) * offset_vector

		line = dwg.add(dwg.polyline(points=[p1, p2, p3], fill="none", stroke='black', stroke_width=1))

		emarker = self.add_e_marker(dwg)
		self.draw_start_arrow(line, emarker)

		dwg.add(dwg.text(textup, insert=text_point_up, fill='black', font_family='sans-serif', font_size=15))
		dwg.add(dwg.text(textdown, insert=text_point_down, fill='black', font_family='sans-serif', font_size=15))

		if element == "weld":
			if self.weld == "Fillet Weld":
				if orientation == "NE":
					self.draw_weld_marker1(dwg, 30, 7.5, line)
				else:
					self.draw_weld_marker2(dwg, 30, 7.5, line)
			else:
				if orientation == "NE":
					self.draw_weld_marker3(dwg, 15, -8.5, line)
				else:
					self.draw_weld_marker4(dwg, 15, 8.5, line)

			if self.stiffener_weld == 1:
				if orientation == "NE":
					self.draw_weld_marker1(dwg, 30, 7.5, line)
				else:
					self.draw_weld_marker2(dwg, 30, 7.5, line)
			else:
				pass

			print "successful"

	def draw_weld_marker1(self, dwg, oriX, oriY, line):
		weldMarker = dwg.marker(insert=(oriX, oriY), size=(15, 15), orient="auto")
		# weldMarker.add(dwg.path(d="M 0 0 L 8 7.5 L 0 15 z", fill='none', stroke='black'))
		weldMarker.add(dwg.path(d="M 15 7.5 L 8 0 L 8 15 z", fill='none', stroke='black'))
		dwg.defs.add(weldMarker)
		self.draw_end_arrow(line, weldMarker)

	def draw_weld_marker2(self, dwg, oriX, oriY, line):
		weldMarker = dwg.marker(insert=(oriX, oriY), size=(15, 15), orient="auto")
		# weldMarker.add(dwg.path(d="M 0 0 L 8 7.5 L 0 15 z", fill='none', stroke='black'))
		weldMarker.add(dwg.path(d="M 0 7.5 L 8 0 L 8 15 z", fill='none', stroke='black'))
		dwg.defs.add(weldMarker)
		self.draw_end_arrow(line, weldMarker)

	def draw_weld_marker3(self, dwg, oriX, oriY, line):
		weldMarker = dwg.marker(insert=(oriX, oriY), size=(15, 15), orient="auto")
		# weldMarker.add(dwg.path(d="M 0 0 L 8 7.5 L 0 15 z", fill='none', stroke='black'))
		weldMarker.add(dwg.path(d="M 0 0 L 0 -7.5 L 7.5 0 ", fill='none', stroke='black'))
		dwg.defs.add(weldMarker)
		self.draw_end_arrow(line, weldMarker)

	def draw_weld_marker4(self, dwg, oriX, oriY, line):
		weldMarker = dwg.marker(insert=(oriX, oriY), size=(15, 15), orient="auto")
		# weldMarker.add(dwg.path(d="M 0 0 L 8 7.5 L 0 15 z", fill='none', stroke='black'))
		weldMarker.add(dwg.path(d="M 0 0 L 0 7.5 L -7.5 0 ", fill='none', stroke='black'))
		dwg.defs.add(weldMarker)
		self.draw_end_arrow(line, weldMarker)

	def save_to_svg(self, filename, view):
		"""

		Args:
			filename: path of the folder
			view: front, top, side views of drawings to be generated

		Returns:
			None

		Note:


		"""
		front_2d = Front_View(self)
		top_2d = Top_View(self)
		side_2d = Side_View(self)
		# extnd_bothway_end_2d_top = ExtendedEnd2DTop(self)
		# extnd_bothway_end_2d_side = ExtendedEnd2DSide(self)

		if view == "Front":
			front_2d.call_Front_View(filename)
			filename = os.path.join(str(self.folder), 'images_html', 'Front.svg')
			front_2d.call_Front_View(filename)
			cairosvg.svg2png(file_obj=filename, write_to=os.path.join(str(self.folder), "images_html", "Front.png"))
		elif view == "Top":
			top_2d.call_Top_View(filename)
			filename = os.path.join(str(self.folder), 'images_html', 'Top.svg')
			top_2d.call_Top_View(filename)
			cairosvg.svg2png(file_obj=filename, write_to=os.path.join(str(self.folder), "images_html", "Top.png"))
		elif view == "Side":
			side_2d.call_Side_View(filename)
			filename = os.path.join(str(self.folder), 'images_html', 'Top.svg')
			side_2d.call_Side_View(filename)
			cairosvg.svg2png(file_obj=filename, write_to=os.path.join(str(self.folder), "images_html", "Side.png"))
		else:
			filename = os.path.join(str(self.folder), 'images_html', 'Front.svg')
			front_2d.call_Front_View(filename)
			cairosvg.svg2png(file_obj=filename, write_to=os.path.join(str(self.folder), "images_html", "Front.png"))

			filename = os.path.join(str(self.folder), 'images_html', 'Top.svg')
			top_2d.call_Top_View(filename)
			cairosvg.svg2png(file_obj=filename, write_to=os.path.join(str(self.folder), "images_html", "Top.png"))

			filename = os.path.join(str(self.folder), 'images_html', 'Side.svg')
			side_2d.call_Side_View(filename)
			cairosvg.svg2png(file_obj=filename, write_to=os.path.join(str(self.folder), "images_html", "Side.png"))


class Front_View(object):
	"""
	Contains functions for generating the front view of the Extended bothway endplate connection.
	"""

	def __init__(self, extnd_common_object):

		self.data_object = extnd_common_object

		# --------------------------------------------------------------------------
		#                               FRONT VIEW (CHANNEL,BEAM AND COLUMNS)
		# --------------------------------------------------------------------------
		# ======================================= Angles =======================================
		if self.data_object.conn_loc == "Leg" or self.data_object.conn_loc == "Back to Back Angles":
			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			ptA2x = ptA1x + self.data_object.member_length
			ptA2y = 0
			self.A2 = np.array([ptA2x, ptA2y])

			if self.data_object.conn_loc == "Leg":
				if self.data_object.weld_oppline > self.data_object.leg_min:
					ptA3x = ptA2x
					ptA3y = ptA2y + self.data_object.leg_max
					self.A3 = np.array([ptA3x, ptA3y])
				else:
					ptA3x = ptA2x
					ptA3y = ptA2y + self.data_object.leg_min
					self.A3 = np.array([ptA3x, ptA3y])
			else:
				if self.data_object.weld_oppline > 2* self.data_object.leg_min:
					ptA3x = ptA2x
					ptA3y = ptA2y + self.data_object.leg_max
					self.A3 = np.array([ptA3x, ptA3y])
				else:
					ptA3x = ptA2x
					ptA3y = ptA2y + self.data_object.leg_min
					self.A3 = np.array([ptA3x, ptA3y])

			ptA4x = ptA3x - self.data_object.member_length
			ptA4y = ptA3y
			self.A4 = np.array([ptA4x, ptA4y])

			ptA5x = ptA4x
			ptA5y = ptA4y - self.data_object.t
			self.A5 = np.array([ptA5x, ptA5y])

			ptA6x = ptA3x
			ptA6y = ptA3y - self.data_object.t
			self.A6 = np.array([ptA6x, ptA6y])

			ptA7x = ptA1x + self.data_object.weld_inline/2
			ptA7y = ptA1y + self.data_object.leg_min/2 - 2 * self.data_object.weld_oppline
			self.A7 = np.array([ptA7x, ptA7y])

			ptA7ax = ptA1x + self.data_object.weld_inline / 2
			ptA7ay = ptA1y
			self.A7a = np.array([ptA7ax, ptA7ay])

			ptA8x = ptA4x + self.data_object.weld_inline / 2
			ptA8y = ptA4y - self.data_object.leg_min / 2 + 2 * self.data_object.weld_oppline
			self.A8 = np.array([ptA8x, ptA8y])

			ptA8ax = ptA4x + self.data_object.weld_inline / 2
			ptA8ay = ptA4y
			self.A8a = np.array([ptA8ax, ptA8ay])

			ptA9x = ptA8x - self.data_object.weld_inline
			ptA9y = ptA8y
			self.A9 = np.array([ptA9x, ptA9y])

			ptA10x = ptA9x
			ptA10y = ptA7y
			self.A10 = np.array([ptA10x, ptA10y])

			if self.data_object.conn_loc == "Back to Back Angles":

				if self.data_object.weld_oppline > 2* self.data_object.leg_min:
					ptA7x = ptA1x + self.data_object.weld_inline / 4
					ptA7y = ptA1y + self.data_object.leg_max / 2 - self.data_object.weld_oppline
					self.A7 = np.array([ptA7x, ptA7y])

					ptA8x = ptA4x + self.data_object.weld_inline / 4
					ptA8y = ptA4y - self.data_object.leg_max / 2 +  self.data_object.weld_oppline
					self.A8 = np.array([ptA8x, ptA8y])
				else:
					ptA7x = ptA1x + self.data_object.weld_inline / 4
					ptA7y = ptA1y + self.data_object.leg_min/ 2 - self.data_object.weld_oppline
					self.A7 = np.array([ptA7x, ptA7y])

					ptA8x = ptA4x + self.data_object.weld_inline / 4
					ptA8y = ptA4y - self.data_object.leg_min / 2 + self.data_object.weld_oppline
					self.A8 = np.array([ptA8x, ptA8y])

				ptA7ax = ptA1x + self.data_object.weld_inline / 4
				ptA7ay = ptA1y
				self.A7a = np.array([ptA7ax, ptA7ay])

				ptA8ax = ptA4x + self.data_object.weld_inline / 4
				ptA8ay = ptA4y
				self.A8a = np.array([ptA8ax, ptA8ay])

				ptA9x = ptA8x - self.data_object.weld_inline/2
				ptA9y = ptA8y
				self.A9 = np.array([ptA9x, ptA9y])

				ptA10x = ptA9x
				ptA10y = ptA7y
				self.A10 = np.array([ptA10x, ptA10y])
			else:
				pass
		# ======================================= Angles =======================================
		elif self.data_object.conn_loc == "Star Angles":

			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			ptA2x = ptA1x + self.data_object.member_length
			ptA2y = 0
			self.A2 = np.array([ptA2x, ptA2y])

			if self.data_object.weld_oppline > 2*self.data_object.leg_min:
				ptA3x = ptA2x
				ptA3y = ptA2y + self.data_object.leg_max
				self.A3 = np.array([ptA3x, ptA3y])
			else:
				ptA3x = ptA2x
				ptA3y = ptA2y + self.data_object.leg_min
				self.A3 = np.array([ptA3x, ptA3y])

			ptA4x = ptA3x - self.data_object.member_length
			ptA4y = ptA3y
			self.A4 = np.array([ptA4x, ptA4y])

			ptA5x = ptA4x
			ptA5y = ptA4y - self.data_object.t
			self.A5 = np.array([ptA5x, ptA5y])

			ptA6x = ptA3x
			ptA6y = ptA3y - self.data_object.t
			self.A6 = np.array([ptA6x, ptA6y])

			ptA7x = ptA4x + self.data_object.weld_inline/4
			ptA7y = ptA4y
			self.A7 = np.array([ptA7x, ptA7y])

			ptA7ax = ptA7x
			ptA7ay = ptA7y - self.data_object.t
			self.A7a = np.array([ptA7ax, ptA7ay])

			ptA8x = ptA7x
			ptA8y = ptA7y + self.data_object.leg_min
			self.A8 = np.array([ptA8x, ptA8y])

			ptA9x = ptA8x - self.data_object.weld_inline/4 + self.data_object.member_length
			ptA9y = ptA8y
			self.A9 = np.array([ptA9x, ptA9y])

			ptA10x = ptA9x - self.data_object.member_length
			ptA10y = ptA9y
			self.A10 = np.array([ptA10x, ptA10y])

			ptA11x = ptA7x
			ptA11y = ptA7y + self.data_object.t
			self.A11 = np.array([ptA11x, ptA11y])

			ptA12x = ptA11x - self.data_object.weld_inline/4 + self.data_object.member_length
			ptA12y = ptA11y
			self.A12 = np.array([ptA12x, ptA12y])

			ptA13x = ptA7x
			ptA13y = ptA7y - self.data_object.weld_oppline
			self.A13 = np.array([ptA13x, ptA13y])

			ptA14x = ptA13x
			ptA14y = ptA13y + 2*self.data_object.weld_oppline
			self.A14 = np.array([ptA14x, ptA14y])

			ptA15x = ptA14x - self.data_object.weld_inline/2
			ptA15y = ptA14y
			self.A15 = np.array([ptA15x, ptA15y])

			ptA16x = ptA15x
			ptA16y = ptA13y
			self.A16 = np.array([ptA16x, ptA16y])

		# ======================================= Channels =======================================
		elif self.data_object.conn_loc == "Back to Back Web" and self.data_object.section_type == "Channels":

			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			ptA2x = ptA1x + (self.data_object.member_length)
			ptA2y = 0
			self.A2 = np.array([ptA2x, ptA2y])

			ptA3x = ptA2x
			ptA3y = ptA2y + self.data_object.member_d
			self.A3 = np.array([ptA3x, ptA3y])

			ptA4x = ptA3x - (self.data_object.member_length)
			ptA4y = ptA3y
			self.A4 = np.array([ptA4x, ptA4y])

			ptA5x = ptA1x
			ptA5y = ptA1y + self.data_object.member_tf
			self.A5 = np.array([ptA5x, ptA5y])

			ptA6x = ptA5x + (self.data_object.member_length)
			ptA6y = ptA5y
			self.A6 = np.array([ptA6x, ptA6y])

			ptA7x = ptA3x
			ptA7y = ptA3y - self.data_object.member_tf
			self.A7 = np.array([ptA7x, ptA7y])

			ptA8x = ptA4x
			ptA8y = ptA4y - self.data_object.member_tf
			self.A8 = np.array([ptA8x, ptA8y])

			ptA9x = ptA1x + self.data_object.weld_inline / 4
			ptA9y = ptA1y - self.data_object.member_d/2
			self.A9 = np.array([ptA9x, ptA9y])

			ptA9ax = ptA1x + self.data_object.weld_inline / 4
			ptA9ay = ptA1y
			self.A9a = np.array([ptA9ax, ptA9ay])

			ptA10x = ptA9x
			ptA10y = ptA9y +  self.data_object.weld_oppline/2 +2* self.data_object.member_d/2
			self.A10 = np.array([ptA10x, ptA10y])

			ptA10ax = ptA4x + self.data_object.weld_inline / 4
			ptA10ay = ptA4y
			self.A10a = np.array([ptA10ax, ptA10ay])

			ptA11x = ptA10x - self.data_object.weld_inline/2
			ptA11y = ptA10y
			self.A11 = np.array([ptA11x, ptA11y])

			ptA12x = ptA9x - self.data_object.weld_inline/2
			ptA12y = ptA9y
			self.A12 = np.array([ptA12x, ptA12y])

		# ======================================= Channels =======================================
		elif self.data_object.conn_loc == "Web" and self.data_object.section_type == "Channels":

			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			ptA2x = ptA1x + (self.data_object.member_length)
			ptA2y = 0
			self.A2 = np.array([ptA2x, ptA2y])

			ptA3x = ptA2x
			ptA3y = ptA2y + self.data_object.member_d
			self.A3 = np.array([ptA3x, ptA3y])

			ptA4x = ptA3x - (self.data_object.member_length)
			ptA4y = ptA3y
			self.A4 = np.array([ptA4x, ptA4y])

			ptA5x = ptA1x
			ptA5y = ptA1y + self.data_object.member_tf
			self.A5 = np.array([ptA5x, ptA5y])

			ptA6x = ptA5x + (self.data_object.member_length)
			ptA6y = ptA5y
			self.A6 = np.array([ptA6x, ptA6y])

			ptA7x = ptA3x
			ptA7y = ptA3y - self.data_object.member_tf
			self.A7 = np.array([ptA7x, ptA7y])

			ptA8x = ptA4x
			ptA8y = ptA4y - self.data_object.member_tf
			self.A8 = np.array([ptA8x, ptA8y])

			ptA9x = ptA1x + self.data_object.weld_inline / 2
			ptA9y = ptA1y + self.data_object.member_d/2- self.data_object.weld_oppline
			self.A9 = np.array([ptA9x, ptA9y])

			ptA9ax = ptA1x + self.data_object.weld_inline / 2
			ptA9ay = ptA1y
			self.A9a = np.array([ptA9ax, ptA9ay])

			ptA10x = ptA9x
			ptA10y = ptA9y + 2* self.data_object.weld_oppline
			self.A10 = np.array([ptA10x, ptA10y])

			ptA10ax = ptA4x + self.data_object.weld_inline / 2
			ptA10ay = ptA4y
			self.A10a = np.array([ptA10ax, ptA10ay])

			ptA11x = ptA10x - self.data_object.weld_inline
			ptA11y = ptA10y
			self.A11 = np.array([ptA11x, ptA11y])

			ptA12x = ptA9x - self.data_object.weld_inline
			ptA12y = ptA9y
			self.A12 = np.array([ptA12x, ptA12y])

		# ======================================= Other than Channels =======================================
		elif self.data_object.conn_loc == "Web" and  self.data_object.section_type != "Channels":
			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			ptA2x = ptA1x + (self.data_object.member_length)
			ptA2y = 0
			self.A2 = np.array([ptA2x, ptA2y])

			ptA3x = ptA2x
			ptA3y = ptA2y + self.data_object.member_d
			self.A3 = np.array([ptA3x, ptA3y])

			ptA4x = ptA3x - (self.data_object.member_length)
			ptA4y = ptA3y
			self.A4 = np.array([ptA4x, ptA4y])

			ptA5x = ptA1x
			ptA5y = ptA1y + self.data_object.member_tf
			self.A5 = np.array([ptA5x, ptA5y])

			ptA6x = ptA5x + (self.data_object.member_length)
			ptA6y = ptA5y
			self.A6 = np.array([ptA6x, ptA6y])

			ptA7x = ptA3x
			ptA7y = ptA3y - self.data_object.member_tf
			self.A7 = np.array([ptA7x, ptA7y])

			ptA8x = ptA4x
			ptA8y = ptA4y - self.data_object.member_tf
			self.A8 = np.array([ptA8x, ptA8y])

			ptA9x = ptA5x
			ptA9y = ptA5y + (self.data_object.member_d- 2*self.data_object.member_tf)/2 - (self.data_object.weld_oppline/2)
			self.A9 = np.array([ptA9x, ptA9y])

			ptA10x = ptA9x + self.data_object.weld_inline/2
			ptA10y = ptA9y
			self.A10 = np.array([ptA10x, ptA10y])

			ptA11x = ptA10x
			ptA11y = ptA10y + self.data_object.weld_oppline
			self.A11 = np.array([ptA11x, ptA11y])

			ptA12x = ptA9x
			ptA12y = ptA11y
			self.A12 = np.array([ptA12x, ptA12y])

			ptA9ax = ptA9x - self.data_object.weld_inline/2
			ptA9ay = ptA9y
			self.A9a = np.array([ptA9ax, ptA9ay])

			ptA12ax = ptA12x - self.data_object.weld_inline/2
			ptA12ay = ptA12y
			self.A12a = np.array([ptA12ax, ptA12ay])

		# ======================================= Except Angles =======================================
		elif self.data_object.conn_loc == "Flange":

			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			ptA2x = ptA1x + (self.data_object.member_length)
			ptA2y = 0
			self.A2 = np.array([ptA2x, ptA2y])

			ptA3x = ptA2x
			ptA3y = ptA2y + self.data_object.member_d
			self.A3 = np.array([ptA3x, ptA3y])

			ptA4x = ptA3x - (self.data_object.member_length)
			ptA4y = ptA3y
			self.A4 = np.array([ptA4x, ptA4y])

			ptA5x = ptA1x
			ptA5y = ptA1y + self.data_object.member_tf
			self.A5 = np.array([ptA5x, ptA5y])

			ptA6x = ptA5x + (self.data_object.member_length)
			ptA6y = ptA5y
			self.A6 = np.array([ptA6x, ptA6y])

			ptA7x = ptA3x
			ptA7y = ptA3y - self.data_object.member_tf
			self.A7 = np.array([ptA7x, ptA7y])

			ptA8x = ptA4x
			ptA8y = ptA4y - self.data_object.member_tf
			self.A8 = np.array([ptA8x, ptA8y])

			ptA9x = ptA1x + self.data_object.weld_inline/4
			ptA9y = ptA1y
			self.A9 = np.array([ptA9x, ptA9y])

			ptA10x = ptA9x
			ptA10y = ptA9y - 12
			self.A10 = np.array([ptA10x, ptA10y])

			ptA11x = ptA10x - self.data_object.weld_inline/2
			ptA11y = ptA10y
			self.A11 = np.array([ptA11x, ptA11y])

			ptA12x = ptA9x - self.data_object.weld_inline / 2
			ptA12y = ptA9y
			self.A12 = np.array([ptA12x, ptA12y])

			ptA13x = ptA4x + self.data_object.weld_inline /4
			ptA13y = ptA4y
			self.A13 = np.array([ptA13x, ptA13y])

			ptA14x = ptA4x + self.data_object.weld_inline /4
			ptA14y = ptA13y + 12
			self.A14 = np.array([ptA14x, ptA14y])

			ptA15x = ptA14x - self.data_object.weld_inline / 2
			ptA15y = ptA14y
			self.A15 = np.array([ptA15x, ptA15y])

			ptA16x = ptA13x - self.data_object.weld_inline / 2
			ptA16y = ptA13y
			self.A16 = np.array([ptA16x, ptA16y])

			# ------------------------------------------  Weld triangle  UP-------------------------------------------
			self.B1 = self.A9
			self.B2 = self.A9 + 12 * np.array([1, 0])
			self.B3 = self.A9 + 12 * np.array([0, -1])

			# ------------------------------------------  Weld triangle  DOWN-------------------------------------------
			self.B4 = self.A13
			self.B5 = self.A13 + 12 * np.array([1, 0])
			self.B6 = self.A13 + 12 * np.array([0, 1])

		else:
			pass

	def call_Front_View(self, filename):
		"""

		Args:
			filename: path of the images to be saved

		Returns:
			Saves the image in the folder

		"""
		# ======================================= Angles =======================================
		if self.data_object.conn_loc == "Leg" or self.data_object.conn_loc =="Back to Back Angles":
			wd = int((self.data_object.member_length) + 1000)
			if self.data_object.weld_oppline > self.data_object.leg_min:
				ht = int(self.data_object.leg_max + 800)
			else:
				ht = int(self.data_object.leg_min + 800)

			dwg = svgwrite.Drawing(filename, size=('100%', '100%'), viewBox=(
				'-600 -400 {} {}').format(wd, ht))
			dwg.add(dwg.polyline(points=[self.A1, self.A2, self.A3, self.A4, self.A1], stroke='blue', fill='none',
								 stroke_width=2.5))
			dwg.add(dwg.line(self.A10, self.A7).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A10, self.A9).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A9, self.A8).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A8, self.A8a).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A7, self.A7a).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A5, self.A6).stroke('blue', width=2.5, linecap='square'))
			pattern = dwg.defs.add(dwg.pattern(id="diagonalHatch", size=(8, 8), patternUnits="userSpaceOnUse",
											   patternTransform="rotate(45 2 2)"))
			pattern.add(dwg.path(d="M 0,1 l 8,0", stroke='#000000', stroke_width=2.5))
			if self.data_object.conn_loc == "Back to Back Angles":
				dwg.add(dwg.rect(insert=(self.A1 - 8 * np.array([0,1])), size=(self.data_object.weld_inline/4, 8), fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))
				dwg.add(dwg.rect(insert=(self.A4), size=(self.data_object.weld_inline/4, 8), fill="url(#diagonalHatch)",
								 stroke='white', stroke_width=1.0))
			else:
				dwg.add(dwg.rect(insert=(self.A1 - 8 * np.array([0, 1])), size=(self.data_object.weld_inline/2, 8),
								 fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))
				dwg.add(dwg.rect(insert=(self.A4), size=(self.data_object.weld_inline/2, 8), fill="url(#diagonalHatch)",
							 stroke='white', stroke_width=1.0))
			dwg.add(dwg.rect(insert=(self.A1- 8 * np.array([1,0])), size=(8, self.data_object.weld_oppline), fill="url(#diagonalHatch)",
							 stroke='white', stroke_width=1.0))

		# ======================================= Angles =======================================
		elif self.data_object.conn_loc == "Star Angles":
			wd = int((self.data_object.member_length) + 1000)
			if self.data_object.weld_oppline > self.data_object.leg_min:
				ht = int(self.data_object.leg_max + 800)
			else:
				ht = int(self.data_object.leg_min + 800)
			dwg = svgwrite.Drawing(filename, size=('100%', '100%'), viewBox=(
				'-600 -400 {} {}').format(wd, ht))
			dwg.add(dwg.polyline(points=[self.A1, self.A2, self.A3, self.A4, self.A1], stroke='blue', fill='none',
								 stroke_width=2.5))

			dwg.add(dwg.line(self.A5, self.A6).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A3,self.A9).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A11, self.A8).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A11, self.A12).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A8, self.A9).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A13, self.A7a).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A8, self.A14).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A13, self.A16).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A15, self.A14).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A15, self.A16).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A4, self.A10).stroke('red', width=2.5, linecap='square').dasharray(dasharray=[5, 5]))
			dwg.add(dwg.line(self.A10, self.A8).stroke('red', width=2.5, linecap='square').dasharray(dasharray=[5, 5]))
			pattern = dwg.defs.add(dwg.pattern(id="diagonalHatch", size=(8, 8), patternUnits="userSpaceOnUse",
											   patternTransform="rotate(45 2 2)"))
			pattern.add(dwg.path(d="M 0,1 l 8,0", stroke='#000000', stroke_width=2.5))
			dwg.add(dwg.rect(insert=(self.A1 - 8 * np.array([0,1])), size=(self.data_object.weld_inline /4, 8), fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))
			dwg.add(dwg.rect(insert=(self.A1 - 8 * np.array([1,0])), size=(8, self.data_object.weld_oppline/2), fill="url(#diagonalHatch)",
							 stroke='white', stroke_width=1.0))
			dwg.add(dwg.rect(insert=self.A4 , size=(self.data_object.weld_inline /4, 8), fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))

		# ======================================= Channels =======================================
		elif self.data_object.conn_loc == "Back to Back Web" or self.data_object.conn_loc == "Web":
			wd = int((self.data_object.member_length) + 1000)
			ht = int(self.data_object.member_d + 800)
			dwg = svgwrite.Drawing(filename, size=('100%', '100%'), viewBox=(
				'-600 -400 {} {}').format(wd, ht))
			dwg.add(dwg.polyline(points=[self.A1, self.A2, self.A3, self.A4, self.A1], stroke='blue', fill='none',
								 stroke_width=2.5))
			dwg.add(dwg.line(self.A5, self.A6).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A8, self.A7).stroke('blue', width=2.5, linecap='square'))
			if self.data_object.section_type == "Channels":
				dwg.add(dwg.line(self.A9, self.A9a).stroke('blue', width=2.5, linecap='square'))

				dwg.add(dwg.line(self.A10a, self.A9a).stroke('red', width=2.5, linecap='square').dasharray(dasharray=[5, 5]))
				dwg.add(dwg.line(self.A10a, self.A10).stroke('blue', width=2.5, linecap='square'))
				dwg.add(dwg.line(self.A10, self.A11).stroke('blue', width=2.5, linecap='square'))
				dwg.add(dwg.line(self.A12, self.A11).stroke('blue', width=2.5, linecap='square'))
				dwg.add(dwg.line(self.A12, self.A9).stroke('blue', width=2.5, linecap='square'))
				pattern = dwg.defs.add(dwg.pattern(id="diagonalHatch", size=(8, 8), patternUnits="userSpaceOnUse",
												   patternTransform="rotate(45 2 2)"))
				pattern.add(dwg.path(d="M 0,1 l 8,0", stroke='#000000', stroke_width=2.5))
				# dwg.add(dwg.rect(insert=(self.A1 - 8 * np.array([0, 1])), size=(self.data_object.weld_inline / 2, 8),
				# 				 fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))
				# dwg.add(dwg.rect(insert=(self.A4), size=(self.data_object.weld_inline / 2, 8), fill="url(#diagonalHatch)",
				# 			 stroke='white', stroke_width=1.0))
				if self.data_object.conn_loc == "Back to Back Web":
					dwg.add(dwg.rect(insert=(self.A1 - 8 * np.array([1, 0])), size=(8, self.data_object.weld_oppline/2),
									 fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))
					dwg.add(dwg.rect(insert=(self.A1 - 8 * np.array([0, 1])), size=(self.data_object.weld_inline / 4, 8),
								 fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))
					dwg.add(dwg.rect(insert=(self.A4), size=(self.data_object.weld_inline / 4, 8),fill="url(#diagonalHatch)",
									 stroke='white', stroke_width=1.0))
				else:
					dwg.add(dwg.rect(insert=(self.A1 - 8 * np.array([1, 0])), size=(8, self.data_object.weld_oppline),
									 fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))
					dwg.add(
						dwg.rect(insert=(self.A1 - 8 * np.array([0, 1])), size=(self.data_object.weld_inline / 2, 8),
								 fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))
					dwg.add(dwg.rect(insert=(self.A4), size=(self.data_object.weld_inline / 2, 8),
									 fill="url(#diagonalHatch)",
									 stroke='white', stroke_width=1.0))
			else:
				dwg.add(dwg.line(self.A9, self.A10).stroke('blue', width=2.5, linecap='square'))
				dwg.add(dwg.line(self.A12, self.A12a).stroke('blue', width=2.5, linecap='square'))
				dwg.add(dwg.line(self.A9, self.A9a).stroke('blue', width=2.5, linecap='square'))
				dwg.add(dwg.line(self.A12a, self.A9a).stroke('blue', width=2.5, linecap='square'))
				dwg.add(dwg.line(self.A10, self.A11).stroke('blue', width=2.5, linecap='square'))
				dwg.add(dwg.line(self.A11, self.A12).stroke('blue', width=2.5, linecap='square'))
				pattern = dwg.defs.add(dwg.pattern(id="diagonalHatch", size=(8, 8), patternUnits="userSpaceOnUse",
												   patternTransform="rotate(45 2 2)"))
				pattern.add(dwg.path(d="M 0,1 l 8,0", stroke='#000000', stroke_width=2.5))
				dwg.add(dwg.rect(insert=(self.A9 - 8 * np.array([0, 1])), size=(self.data_object.weld_inline / 2, 8),
								 fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))
				dwg.add(dwg.rect(insert=(self.A12), size=(self.data_object.weld_inline / 2, 8), fill="url(#diagonalHatch)",
							 stroke='white', stroke_width=1.0))
				dwg.add(dwg.rect(insert=(self.A10), size=(8, self.data_object.weld_oppline),
								 fill="url(#diagonalHatch)",
								 stroke='white', stroke_width=1.0))

		# ======================================= Other than Angles =======================================
		elif self.data_object.conn_loc == "Flange":
			wd = int((self.data_object.member_length) + 1000)
			ht = int(self.data_object.member_d + 800)
			dwg = svgwrite.Drawing(filename, size=('100%', '100%'), viewBox=(
				'-600 -400 {} {}').format(wd, ht))
			dwg.add(dwg.polyline(points=[self.A1, self.A2, self.A3, self.A4, self.A1], stroke='blue', fill='none', stroke_width=2.5))
			dwg.add(dwg.line(self.A5, self.A6).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A8, self.A7).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A9, self.A10).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A10, self.A11).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A11, self.A12).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A1, self.A12).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A13, self.A14).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A14, self.A15).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A15, self.A16).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A16, self.A4).stroke('blue', width=2.5, linecap='square'))
			pattern = dwg.defs.add(dwg.pattern(id="diagonalHatch", size=(8, 8), patternUnits="userSpaceOnUse",
											   patternTransform="rotate(45 2 2)"))
			pattern.add(dwg.path(d="M 0,1 l 8,0", stroke='#000000', stroke_width=2.5))
			dwg.add(dwg.rect(insert=(self.A1 - 8 * np.array([0, 1])), size=(self.data_object.weld_inline / 4, 8),
							 fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))
			dwg.add(dwg.rect(insert=(self.A4), size=(self.data_object.weld_inline / 4, 8), fill="url(#diagonalHatch)",
						 stroke='white', stroke_width=1.0))
			dwg.add(dwg.polyline(points=[self.B1, self.B2, self.B3, self.B1], stroke='black', fill='red',
								 stroke_width=2.5))
			dwg.add(dwg.polyline(points=[self.B4, self.B5, self.B6, self.B4], stroke='black', fill='red',
								 stroke_width=2.5))

		else:
			pass
		# ------------------------------------------  Beam Designation Labeling  -------------------------------------------
		point = self.A1
		theta = 30
		offset = 100
		textup = str(self.data_object.beam_designation)
		textdown = " "
		element = " "
		self.data_object.draw_oriented_arrow(dwg, point, theta, "NW", offset, textup, textdown, element)

		# ------------------------------------------  Beam Designation Labeling  -------------------------------------------
		pt_a1 = self.A1 + (300) * np.array([0, -1])- (200) * np.array([1, 0])
		pt_b1 = pt_a1 + (50 * np.array([0, 1]))
		txt_1 = pt_b1 + (20 * np.array([-1, 0])) + (75 * np.array([0, 1]))
		text = "A"
		self.data_object.draw_cross_section(dwg, pt_a1, pt_b1, txt_1, text)

		pt_a2 = self.A2 + (300) * np.array([0, -1]) + (200) * np.array([1, 0])
		pt_b2 = pt_a2 + (50 * np.array([0, 1]))
		txt_2 = pt_b2 + (20 * np.array([-1, 0])) + (75 * np.array([0, 1]))
		self.data_object.draw_cross_section(dwg, pt_a2, pt_b2, txt_2, text)

		dwg.add(dwg.line(pt_a1, pt_a2).stroke('black', width=1.5, linecap='square'))

		# ------------------------------------------  Dimensions Labeling  -------------------------------------------------

		if self.data_object.conn_loc == "Star Angles":
			ptx1 = self.A4
			pty1 = ptx1 + (150) * np.array([0, 1])
			self.data_object.draw_faint_line(ptx1, pty1, dwg)

			ptx2 = self.A3
			pty2 = ptx2 + (150) * np.array([0, 1])
			self.data_object.draw_faint_line(ptx2, pty2, dwg)

			point1 = ptx2 - (self.data_object.member_length) * np.array([1, 0])
			params = {"offset": -150, "textoffset": 5, "lineori": "right", "endlinedim": 10, "arrowlen": 20}
			self.data_object.draw_dimension_outer_arrow(dwg, ptx2, point1, str(self.data_object.act_member_length),
														params)
		else:
			ptx1 = self.A4
			pty1 = ptx1 + (100) * np.array([0, 1])
			self.data_object.draw_faint_line(ptx1, pty1, dwg)

			ptx2 = self.A3
			pty2 = ptx2 + (100) * np.array([0, 1])
			self.data_object.draw_faint_line(ptx2, pty2, dwg)

			point1 = ptx2 - (self.data_object.member_length) * np.array([1, 0])
			params = {"offset": -100, "textoffset": 5, "lineori": "right", "endlinedim": 10, "arrowlen": 20}
			self.data_object.draw_dimension_outer_arrow(dwg, ptx2, point1, str(self.data_object.act_member_length),
														params)



		# ptx3 = self.A4
		# pty1 = ptx1 + (100) * np.array([0, 1])
		# self.data_object.draw_faint_line(ptx1, pty1, dwg)
		#
		# ptx2 = self.A3
		# pty2 = ptx2 + (100) * np.array([0, 1])
		# self.data_object.draw_faint_line(ptx2, pty2, dwg)
		#
		# point1 = ptx2 - (self.data_object.member_length) * np.array([1, 0])
		# params = {"offset": -150, "textoffset": 10, "lineori": "right", "endlinedim": 10, "arrowlen": 20}
		# self.data_object.draw_dimension_outer_arrow(dwg, ptx2, point1, str(self.data_object.act_member_length),
		# 											params)

		# ------------------------------------------  View details-------------------------------------------
		ptx = self.A4 + ((self.data_object.member_length/2)-300)* np.array([1, 0]) + 250 * np.array([0, 1])
		dwg.add(dwg.text('Front view (Sec C-C) ', insert=ptx, fill='black', font_family="sans-serif", font_size=35))
		ptx1 = ptx + 40 * np.array([0, 1])
		dwg.add(dwg.text('(All dimensions are in "mm")', insert=ptx1, fill='black', font_family="sans-serif", font_size=35))
		dwg.save()


class Top_View(object):
	"""
	Contains functions for generating the top view of the Extended bothway endplate connection.

	"""

	def __init__(self, extnd_common_object):
		self.data_object = extnd_common_object
		# -------------------------------------------------------------------------------------------------
		#                                           TOP VIEW
		# -------------------------------------------------------------------------------------------------
		# ============================================ Angles =============================================
		if self.data_object.conn_loc == "Leg":
			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			ptA2x = ptA1x + self.data_object.member_length
			ptA2y = 0
			self.A2 = np.array([ptA2x, ptA2y])

			if self.data_object.weld_oppline > self.data_object.leg_min:
				ptA3x = ptA2x
				ptA3y = ptA2y + self.data_object.leg_max
				self.A3 = np.array([ptA3x, ptA3y])
			else:
				ptA3x = ptA2x
				ptA3y = ptA2y + self.data_object.leg_min
				self.A3 = np.array([ptA3x, ptA3y])

			ptA4x = ptA1x
			ptA4y = ptA3y
			self.A4 = np.array([ptA4x, ptA4y])

			ptA5x = ptA1x
			ptA5y = ptA1y + self.data_object.t
			self.A5 = np.array([ptA5x, ptA5y])

			ptA6x = ptA2x
			ptA6y = ptA2y + self.data_object.t
			self.A6 = np.array([ptA6x, ptA6y])

			ptA7x = ptA1x + self.data_object.weld_inline/2
			ptA7y = ptA1y
			self.A7 = np.array([ptA7x, ptA7y])

			ptA8x = ptA7x
			ptA8y = ptA7y - 12
			self.A8 = np.array([ptA8x, ptA8y])

			ptA9x = ptA8x - self.data_object.weld_inline
			ptA9y = ptA8y
			self.A9= np.array([ptA9x, ptA9y])

			ptA10x = ptA9x
			ptA10y = ptA9y + 12
			self.A10 = np.array([ptA10x, ptA10y])

		elif self.data_object.conn_loc == "Back to Back Angles":

			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			ptA2x = ptA1x + self.data_object.member_length
			ptA2y = 0
			self.A2 = np.array([ptA2x, ptA2y])

			if self.data_object.weld_oppline > 2*self.data_object.leg_min:
				ptA3x = ptA2x
				ptA3y = ptA2y + self.data_object.leg_max
				self.A3 = np.array([ptA3x, ptA3y])
			else:
				ptA3x = ptA2x
				ptA3y = ptA2y + self.data_object.leg_min
				self.A3 = np.array([ptA3x, ptA3y])

			ptA4x = ptA1x
			ptA4y = ptA3y
			self.A4 = np.array([ptA4x, ptA4y])

			ptA5x = ptA1x
			ptA5y = ptA1y + self.data_object.t
			self.A5 = np.array([ptA5x, ptA5y])

			ptA6x = ptA2x
			ptA6y = ptA2y + self.data_object.t
			self.A6 = np.array([ptA6x, ptA6y])

			ptA7x = ptA1x + self.data_object.weld_inline / 4
			ptA7y = ptA1y
			self.A7 = np.array([ptA7x, ptA7y])

			ptA8x = ptA7x
			ptA8y = ptA7y - 12
			self.A8 = np.array([ptA8x, ptA8y])

			ptA9x = ptA8x - self.data_object.weld_inline/2
			ptA9y = ptA8y
			self.A9 = np.array([ptA9x, ptA9y])

			ptA10x = ptA9x
			ptA10y = ptA9y + 12
			self.A10 = np.array([ptA10x, ptA10y])

			ptA11x = ptA1x
			ptA11y = ptA1y - 12
			self.A11 = np.array([ptA11x, ptA11y])

			ptA12x = ptA11x + self.data_object.member_length
			ptA12y = ptA11y
			self.A12 = np.array([ptA12x, ptA12y])

			if self.data_object.weld_oppline > 2*self.data_object.leg_min:
				ptA13x = ptA12x
				ptA13y = ptA12y - self.data_object.leg_max
				self.A13 = np.array([ptA13x, ptA13y])
			else:
				ptA13x = ptA12x
				ptA13y = ptA12y - self.data_object.leg_min
				self.A13 = np.array([ptA13x, ptA13y])

			ptA14x = ptA11x
			ptA14y = ptA13y
			self.A14 = np.array([ptA14x, ptA14y])

			ptA15x = ptA11x
			ptA15y = ptA11y - self.data_object.t
			self.A15 = np.array([ptA15x, ptA15y])

			ptA16x = ptA12x
			ptA16y = ptA12y - self.data_object.t
			self.A16 = np.array([ptA16x, ptA16y])

		elif self.data_object.conn_loc == "Star Angles":

			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			ptA2x = ptA1x + self.data_object.member_length
			ptA2y = 0
			self.A2 = np.array([ptA2x, ptA2y])

			if self.data_object.weld_oppline > 2* self.data_object.leg_min:
				ptA3x = ptA2x
				ptA3y = ptA2y + self.data_object.leg_max
				self.A3 = np.array([ptA3x, ptA3y])
			else:
				ptA3x = ptA2x
				ptA3y = ptA2y + self.data_object.leg_min
				self.A3 = np.array([ptA3x, ptA3y])

			ptA4x = ptA1x
			ptA4y = ptA3y
			self.A4 = np.array([ptA4x, ptA4y])

			ptA5x = ptA1x
			ptA5y = ptA1y + self.data_object.t
			self.A5 = np.array([ptA5x, ptA5y])

			ptA6x = ptA2x
			ptA6y = ptA2y + self.data_object.t
			self.A6 = np.array([ptA6x, ptA6y])

			ptA7x = ptA1x + self.data_object.weld_inline / 4
			ptA7y = ptA1y
			self.A7 = np.array([ptA7x, ptA7y])

			ptA8x = ptA7x
			ptA8y = ptA7y - 12
			self.A8 = np.array([ptA8x, ptA8y])

			ptA9x = ptA8x - self.data_object.weld_inline/2
			ptA9y = ptA8y
			self.A9 = np.array([ptA9x, ptA9y])

			ptA10x = ptA9x
			ptA10y = ptA9y + 12
			self.A10 = np.array([ptA10x, ptA10y])

			ptA11x = ptA1x
			ptA11y = ptA1y - 12
			self.A11 = np.array([ptA11x, ptA11y])

			ptA12x = ptA11x + self.data_object.member_length
			ptA12y = ptA11y
			self.A12 = np.array([ptA12x, ptA12y])

			if self.data_object.weld_oppline > 2* self.data_object.leg_min:
				ptA13x = ptA12x
				ptA13y = ptA12y - self.data_object.leg_max
				self.A13 = np.array([ptA13x, ptA13y])
			else:
				ptA13x = ptA12x
				ptA13y = ptA12y - self.data_object.leg_min
				self.A13 = np.array([ptA13x, ptA13y])

			ptA14x = ptA11x
			ptA14y = ptA13y
			self.A14 = np.array([ptA14x, ptA14y])

			ptA15x = ptA11x
			ptA15y = ptA11y - self.data_object.t
			self.A15 = np.array([ptA15x, ptA15y])

			ptA16x = ptA12x
			ptA16y = ptA12y - self.data_object.t
			self.A16 = np.array([ptA16x, ptA16y])

		elif (self.data_object.conn_loc == "Web" and (self.data_object.section_type == "Columns" or self.data_object.section_type == "Beams")):

			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			ptA2x = ptA1x + self.data_object.member_length
			ptA2y = 0
			self.A2 = np.array([ptA2x, ptA2y])

			# if self.data_object.weld_oppline > self.data_object.leg_min:
			ptA3x = ptA2x
			ptA3y = ptA2y + self.data_object.member_B
			self.A3 = np.array([ptA3x, ptA3y])
			# else:
			# 	ptA3x = ptA2x
			# 	ptA3y = ptA2y + self.data_object.leg_min
			# 	self.A3 = np.array([ptA3x, ptA3y])

			ptA4x = ptA1x
			ptA4y = ptA3y
			self.A4 = np.array([ptA4x, ptA4y])

			ptA5x = ptA1x
			ptA5y = ptA1y + self.data_object.member_B/2 - self.data_object.member_tw/2
			self.A5 = np.array([ptA5x, ptA5y])

			ptA6x = ptA5x + self.data_object.member_length
			ptA6y = ptA5y
			self.A6 = np.array([ptA6x, ptA6y])

			ptA7x = ptA6x
			ptA7y = ptA5y + self.data_object.member_tw
			self.A7 = np.array([ptA7x, ptA7y])

			ptA8x = ptA5x
			ptA8y = ptA7y
			self.A8 = np.array([ptA8x, ptA8y])

			ptA9x = ptA8x + self.data_object.weld_inline/2
			ptA9y = ptA8y
			self.A9 = np.array([ptA9x, ptA9y])

			ptA10x = ptA9x
			ptA10y = ptA9y + 12
			self.A10 = np.array([ptA10x, ptA10y])

			ptA11x = ptA8x
			ptA11y = ptA8y + 12
			self.A11 = np.array([ptA11x, ptA11y])

			ptA12x = ptA11x - self.data_object.weld_inline/2
			ptA12y = ptA11y
			self.A12 = np.array([ptA12x, ptA12y])

			# if self.data_object.weld_oppline > self.data_object.leg_min:
			ptA13x = ptA12x
			ptA13y = ptA12y - 12
			self.A13 = np.array([ptA13x, ptA13y])

		elif self.data_object.conn_loc == "Web" or self.data_object.conn_loc == "Back to Back Web" and self.data_object.section_type == "Channels":

			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			ptA2x = ptA1x + self.data_object.member_length
			ptA2y = 0
			self.A2 = np.array([ptA2x, ptA2y])

			ptA3x = ptA2x
			ptA3y = ptA2y + self.data_object.member_B
			self.A3 = np.array([ptA3x, ptA3y])

			ptA4x = ptA1x
			ptA4y = ptA3y
			self.A4 = np.array([ptA4x, ptA4y])

			ptA5x = ptA1x
			ptA5y = ptA1y + self.data_object.member_tw
			self.A5 = np.array([ptA5x, ptA5y])

			ptA6x = ptA5x + self.data_object.member_length
			ptA6y = ptA5y
			self.A6 = np.array([ptA6x, ptA6y])

			if self.data_object.conn_loc == "Back to Back Web":
				ptA9x = ptA1x + self.data_object.weld_inline/4
				ptA9y = ptA1y
				self.A9 = np.array([ptA9x, ptA9y])

				ptA10x = ptA9x
				ptA10y = ptA9y - 12
				self.A10 = np.array([ptA10x, ptA10y])

				ptA11x = ptA10x - self.data_object.weld_inline/2
				ptA11y = ptA10y
				self.A11 = np.array([ptA11x, ptA11y])

				ptA12x = ptA11x
				ptA12y = ptA11y + 12
				self.A12 = np.array([ptA12x, ptA12y])

				ptA13x = ptA12x
				ptA13y = ptA12y - 12
				self.A13 = np.array([ptA13x, ptA13y])
			else:
				ptA9x = ptA1x + self.data_object.weld_inline / 2
				ptA9y = ptA1y
				self.A9 = np.array([ptA9x, ptA9y])

				ptA10x = ptA9x
				ptA10y = ptA9y - 12
				self.A10 = np.array([ptA10x, ptA10y])

				ptA11x = ptA10x - self.data_object.weld_inline
				ptA11y = ptA10y
				self.A11 = np.array([ptA11x, ptA11y])

				ptA12x = ptA11x
				ptA12y = ptA11y + 12
				self.A12 = np.array([ptA12x, ptA12y])

				ptA13x = ptA12x
				ptA13y = ptA12y - 12
				self.A13 = np.array([ptA13x, ptA13y])


			if self.data_object.conn_loc == "Back to Back Web":
				ptA14x = ptA1x
				ptA14y = ptA1y - 12
				self.A14 = np.array([ptA14x, ptA14y])

				ptA15x = ptA1x + self.data_object.member_length
				ptA15y = ptA14y
				self.A15 = np.array([ptA15x, ptA15y])

				ptA16x = ptA15x
				ptA16y = ptA15y - self.data_object.member_B
				self.A16 = np.array([ptA16x, ptA16y])

				ptA17x = ptA14x
				ptA17y = ptA16y
				self.A17 = np.array([ptA17x, ptA17y])

				ptA18x = ptA14x
				ptA18y = ptA14y - self.data_object.member_tw
				self.A18 = np.array([ptA18x, ptA18y])

				ptA19x = ptA15x
				ptA19y = ptA15y - self.data_object.member_tw
				self.A19 = np.array([ptA19x, ptA19y])

		elif self.data_object.conn_loc == "Flange" and (self.data_object.section_type == "Columns" or self.data_object.section_type == "Beams"):

			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			ptA2x = ptA1x + self.data_object.member_length
			ptA2y = 0
			self.A2 = np.array([ptA2x, ptA2y])

			ptA3x = ptA2x
			ptA3y = ptA2y + self.data_object.member_B
			self.A3 = np.array([ptA3x, ptA3y])

			ptA4x = ptA1x
			ptA4y = ptA3y
			self.A4 = np.array([ptA4x, ptA4y])

			ptA5x = ptA1x
			ptA5y = ptA1y + self.data_object.member_B / 2 - self.data_object.member_tw / 2
			self.A5 = np.array([ptA5x, ptA5y])

			ptA6x = ptA5x + self.data_object.member_length
			ptA6y = ptA5y
			self.A6 = np.array([ptA6x, ptA6y])

			ptA7x = ptA6x
			ptA7y = ptA5y + self.data_object.member_tw
			self.A7 = np.array([ptA7x, ptA7y])

			ptA8x = ptA5x
			ptA8y = ptA7y
			self.A8 = np.array([ptA8x, ptA8y])

			ptA9x = ptA1x
			ptA9y = ptA1y - self.data_object.weld_oppline/4 + self.data_object.member_B/2
			self.A9 = np.array([ptA9x, ptA9y])

			ptA10x = ptA9x + self.data_object.weld_inline/4
			ptA10y = ptA9y
			self.A10 = np.array([ptA10x, ptA10y])

			ptA11x = ptA10x
			ptA11y = ptA10y + self.data_object.weld_oppline/2
			self.A11 = np.array([ptA11x, ptA11y])

			ptA12x = ptA11x - self.data_object.weld_inline/4
			ptA12y = ptA11y
			self.A12 = np.array([ptA12x, ptA12y])

			ptA13x = ptA12x - self.data_object.weld_inline / 4
			ptA13y = ptA12y
			self.A13 = np.array([ptA13x, ptA13y])

			ptA14x = ptA13x
			ptA14y = ptA13y - self.data_object.weld_oppline/2
			self.A14 = np.array([ptA14x, ptA14y])

		elif self.data_object.conn_loc == "Flange" and self.data_object.section_type == "Channels":

			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			ptA2x = ptA1x + self.data_object.member_length
			ptA2y = 0
			self.A2 = np.array([ptA2x, ptA2y])

			ptA3x = ptA2x
			ptA3y = ptA2y + self.data_object.member_B
			self.A3 = np.array([ptA3x, ptA3y])

			ptA4x = ptA1x
			ptA4y = ptA3y
			self.A4 = np.array([ptA4x, ptA4y])

			ptA5x = ptA1x
			ptA5y = ptA1y + self.data_object.member_tw
			self.A5 = np.array([ptA5x, ptA5y])

			ptA6x = ptA5x + self.data_object.member_length
			ptA6y = ptA5y
			self.A6 = np.array([ptA6x, ptA6y])

			ptA9x = ptA1x
			ptA9y = ptA1y - self.data_object.weld_oppline / 4 + self.data_object.member_B / 2
			self.A9 = np.array([ptA9x, ptA9y])

			ptA10x = ptA9x + self.data_object.weld_inline / 4
			ptA10y = ptA9y
			self.A10 = np.array([ptA10x, ptA10y])

			ptA11x = ptA10x
			ptA11y = ptA10y + self.data_object.weld_oppline / 2
			self.A11 = np.array([ptA11x, ptA11y])

			ptA12x = ptA11x - self.data_object.weld_inline / 4
			ptA12y = ptA11y
			self.A12 = np.array([ptA12x, ptA12y])

			ptA13x = ptA12x - self.data_object.weld_inline / 4
			ptA13y = ptA12y
			self.A13 = np.array([ptA13x, ptA13y])

			ptA14x = ptA13x
			ptA14y = ptA13y - self.data_object.weld_oppline / 2
			self.A14 = np.array([ptA14x, ptA14y])

	def call_Top_View(self, filename):
		"""

		Args:
			filename: path of the images to be saved

		Returns:
			Saves the image in the folder

		"""
		if self.data_object.conn_loc == "Leg":
			wd = int((self.data_object.member_length) + 1000)
			if self.data_object.weld_oppline > self.data_object.leg_min:
				ht = int(self.data_object.leg_max + 800)
			else:
				ht = int(self.data_object.leg_min + 800)
			dwg = svgwrite.Drawing(filename, size=('100%', '100%'), viewBox=(
				'-600 -400 {} {}').format(wd, ht))
			dwg.add(dwg.polyline(points=[self.A1, self.A2, self.A3, self.A4, self.A1], stroke='blue', fill='none',
								 stroke_width=2.5))
			dwg.add(dwg.line(self.A5, self.A6).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A8, self.A7).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A8, self.A9).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A10, self.A9).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A10, self.A1).stroke('blue', width=2.5, linecap='square'))
			pattern = dwg.defs.add(dwg.pattern(id="diagonalHatch", size=(8, 8), patternUnits="userSpaceOnUse",
											   patternTransform="rotate(45 2 2)"))
			pattern.add(dwg.path(d="M 0,1 l 8,0", stroke='#000000', stroke_width=2.5))
			dwg.add(dwg.rect(insert=(self.A1), size=(self.data_object.weld_inline / 2, self.data_object.t ),
							 fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))


		elif self.data_object.conn_loc == "Back to Back Angles":
			wd = int((self.data_object.member_length) + 1000)
			if self.data_object.weld_oppline > self.data_object.leg_min:
				ht = int(self.data_object.leg_max + 800)
			else:
				ht = int(self.data_object.leg_min + 800)
			dwg = svgwrite.Drawing(filename, size=('100%', '100%'), viewBox=(
				'-600 -400 {} {}').format(wd, ht))
			dwg.add(dwg.polyline(points=[self.A1, self.A2, self.A3, self.A4, self.A1], stroke='blue', fill='none',
								 stroke_width=2.5))
			dwg.add(dwg.line(self.A5, self.A6).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A8, self.A7).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A8, self.A9).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A10, self.A9).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A10, self.A1).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A11, self.A14).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A14, self.A13).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A12, self.A13).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A12, self.A8).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A15, self.A16).stroke('blue', width=2.5, linecap='square'))
			pattern = dwg.defs.add(dwg.pattern(id="diagonalHatch", size=(8, 8), patternUnits="userSpaceOnUse",
											   patternTransform="rotate(45 2 2)"))
			pattern.add(dwg.path(d="M 0,1 l 8,0", stroke='#000000', stroke_width=2.5))
			dwg.add(dwg.rect(insert=(self.A1), size=(self.data_object.weld_inline /4, self.data_object.t),
							 fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))
			dwg.add(dwg.rect(insert=(self.A11 -(self.data_object.t* np.array([0, 1]))), size=(self.data_object.weld_inline / 4, self.data_object.t),
							 fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))

		elif self.data_object.conn_loc == "Star Angles":
			wd = int((self.data_object.member_length) + 1000)
			if self.data_object.weld_oppline > self.data_object.leg_min:
				ht = int(self.data_object.leg_max + 800)
			else:
				ht = int(self.data_object.leg_min + 800)
			dwg = svgwrite.Drawing(filename, size=('100%', '100%'), viewBox=(
				'-600 -400 {} {}').format(wd, ht))
			dwg.add(dwg.polyline(points=[self.A1, self.A2, self.A3, self.A4, self.A1], stroke='blue', fill='none',
								 stroke_width=2.5))
			dwg.add(dwg.line(self.A5, self.A6).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A8, self.A7).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A8, self.A9).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A10, self.A9).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A10, self.A1).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A11, self.A14).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A14, self.A13).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A12, self.A13).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A12, self.A8).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A15, self.A16).stroke('blue', width=2.5, linecap='square').dasharray(dasharray=[5, 5]))
			pattern = dwg.defs.add(dwg.pattern(id="diagonalHatch", size=(8, 8), patternUnits="userSpaceOnUse",
											   patternTransform="rotate(45 2 2)"))
			pattern.add(dwg.path(d="M 0,1 l 8,0", stroke='#000000', stroke_width=2.5))
			dwg.add(dwg.rect(insert=(self.A1), size=(self.data_object.weld_inline / 4, self.data_object.t),
							 fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))
			dwg.add(dwg.rect(insert=(self.A11 - (self.data_object.t * np.array([0, 1]))),
							 size=(self.data_object.weld_inline / 4, self.data_object.t),
							 fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))

		elif self.data_object.conn_loc == "Web" or self.data_object.conn_loc == "Back to Back Web":
			if self.data_object.conn_loc == "Back to Back Web":
				wd = int((self.data_object.member_length) + 750)
				ht = int(self.data_object.member_d + 500)
				dwg = svgwrite.Drawing(filename, size=('100%', '100%'), viewBox=(
					'-375 -300 {} {}').format(wd, ht))
			else:
				wd = int((self.data_object.member_length) + 1000)
				ht = int(self.data_object.member_d + 800)
				dwg = svgwrite.Drawing(filename, size=('100%', '100%'), viewBox=(
					'-600 -400 {} {}').format(wd, ht))

			dwg.add(dwg.polyline(points=[self.A1, self.A2, self.A3, self.A4, self.A1], stroke='blue', fill='none',
								 stroke_width=2.5))

			dwg.add(dwg.line(self.A5, self.A6).stroke('blue', width=2.5, linecap='square').dasharray(dasharray=[5, 5]))
			if self.data_object.section_type == "Columns" or self.data_object.section_type == "Beams":
				dwg.add(dwg.line(self.A7, self.A8).stroke('blue', width=2.5, linecap='square').dasharray(dasharray=[5, 5]))
				dwg.add(dwg.line(self.A13, self.A8).stroke('blue', width=2.5, linecap='square'))
				dwg.add(dwg.line(self.A9, self.A10).stroke('blue', width=2.5, linecap='square').dasharray(dasharray=[5, 5]))
				dwg.add(dwg.line(self.A10, self.A11).stroke('blue', width=2.5, linecap='square').dasharray(dasharray=[5, 5]))
				dwg.add(dwg.line(self.A12, self.A13).stroke('blue', width=2.5, linecap='square'))
				dwg.add(dwg.line(self.A11, self.A12).stroke('blue', width=2.5, linecap='square'))
			else:
				# dwg.add(dwg.line(self.A13, self.A5).stroke('blue', width=2.5, linecap='square'))
				dwg.add(dwg.line(self.A9, self.A10).stroke('blue', width=2.5, linecap='square'))
				dwg.add(dwg.line(self.A10, self.A11).stroke('blue', width=2.5, linecap='square'))
				dwg.add(dwg.line(self.A11, self.A12).stroke('blue', width=2.5, linecap='square'))
				dwg.add(dwg.line(self.A1, self.A12).stroke('blue', width=2.5, linecap='square'))
				if self.data_object.conn_loc == "Back to Back Web":
					dwg.add(dwg.line(self.A14, self.A15).stroke('blue', width=2.5, linecap='square'))
					dwg.add(dwg.line(self.A15, self.A16).stroke('blue', width=2.5, linecap='square'))
					dwg.add(dwg.line(self.A16, self.A17).stroke('blue', width=2.5, linecap='square'))
					dwg.add(dwg.line(self.A17, self.A14).stroke('blue', width=2.5, linecap='square'))
					dwg.add(dwg.line(self.A18, self.A19).stroke('blue', width=2.5, linecap='square').dasharray(
						dasharray=[5, 5]))
				else:
					pass
		elif self.data_object.conn_loc == "Flange":
			wd = int((self.data_object.member_length) + 1000)
			ht = int(self.data_object.member_d + 800)
			dwg = svgwrite.Drawing(filename, size=('100%', '100%'), viewBox=(
				'-600 -400 {} {}').format(wd, ht))
			# dwg.add(dwg.polyline(points=[self.A1, self.A2, self.A3, self.A4, self.A1], stroke='blue', fill='none',
			# 					 stroke_width=2.5))
			dwg.add(dwg.line(self.A1, self.A2).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A2, self.A3).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A3, self.A4).stroke('blue', width=2.5, linecap='square'))
			# dwg.add(dwg.line(self.A1, self.A7).stroke('blue', width=2.5, linecap='square'))
			# dwg.add(dwg.line(self.A4, self.A10).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A5, self.A6).stroke('blue', width=2.5, linecap='square').dasharray(dasharray=[5, 5]))
			if self.data_object.section_type != "Channels":
				dwg.add(dwg.line(self.A7, self.A8).stroke('blue', width=2.5, linecap='square').dasharray(dasharray=[5, 5]))
			else:
				pass
			dwg.add(dwg.line(self.A9, self.A1).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A4, self.A12).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A9, self.A12).stroke('blue', width=2.5, linecap='square').dasharray(dasharray=[5, 5]))
			dwg.add(dwg.line(self.A9, self.A10).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A10, self.A11).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A11, self.A12).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A12, self.A13).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A13, self.A14).stroke('blue', width=2.5, linecap='square'))
			dwg.add(dwg.line(self.A14, self.A9).stroke('blue', width=2.5, linecap='square'))
			pattern = dwg.defs.add(dwg.pattern(id="diagonalHatch", size=(8, 8), patternUnits="userSpaceOnUse",
											   patternTransform="rotate(45 2 2)"))
			pattern.add(dwg.path(d="M 0,1 l 8,0", stroke='#000000', stroke_width=2.5))
			dwg.add(dwg.rect(insert=(self.A9 - (12 * np.array([0, 1]))), size=(self.data_object.weld_inline / 4, 8),
							 fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))
			dwg.add(dwg.rect(insert=(self.A12),
							 size=(self.data_object.weld_inline / 4, 8),
							 fill="url(#diagonalHatch)", stroke='white', stroke_width=1.0))
			dwg.add(dwg.rect(insert=(self.A10), size=(8, self.data_object.weld_oppline/2),
							 fill="url(#diagonalHatch)",
							 stroke='white', stroke_width=1.0))

		# ------------------------------------------  Primary Beam 1& 2 -------------------------------------------
		point = self.A1
		theta = 60
		offset = 100
		textdown = " "
		textup =  str(self.data_object.beam_designation)
		element = " "
		self.data_object.draw_oriented_arrow(dwg, point, theta, "NW", offset, textup, textdown, element)

		# ------------------------------------------Dimension Labelling -------------------------------------------
		if self.data_object.section_type == "Angles":
			ptx1 = self.A1
			pty1 = ptx1 + (200) * np.array([0, -1])
			self.data_object.draw_faint_line(ptx1, pty1, dwg)

			ptx2 = self.A7
			pty2 = ptx2 + (200) * np.array([0, -1])
			self.data_object.draw_faint_line(ptx2, pty2, dwg)

			if self.data_object.conn_loc == "Leg":
				point1 = ptx2 - (self.data_object.weld_inline/2) * np.array([1, 0])
				params = {"offset": -200, "textoffset": 10, "lineori": "left", "endlinedim": 10, "arrowlen": 20}
				self.data_object.draw_dimension_outer_arrow(dwg,ptx2, point1,str(self.data_object.weld_inline/2),params)
			elif self.data_object.conn_loc == "Back to Back Angles":
				point1 = ptx2 - (self.data_object.weld_inline / 4) * np.array([1, 0])
				params = {"offset": -200, "textoffset": 10, "lineori": "left", "endlinedim": 10, "arrowlen": 20}
				self.data_object.draw_dimension_outer_arrow(dwg, ptx2, point1, str(self.data_object.weld_inline/4),params)
			elif self.data_object.conn_loc == "Star Angles":
				point1 = ptx2 - (self.data_object.weld_inline / 4) * np.array([1, 0])
				params = {"offset": -200, "textoffset": 10, "lineori": "left", "endlinedim": 10, "arrowlen": 20}
				self.data_object.draw_dimension_outer_arrow(dwg, ptx2, point1, str(self.data_object.weld_inline/4),params)

		elif self.data_object.section_type == "Beams" or self.data_object.section_type == "Columns"or self.data_object.section_type == "Channels" :
			if self.data_object.conn_loc != "Web":
				ptx1 = self.A9
				pty1 = ptx1 + (100) * np.array([0, -1])
				self.data_object.draw_faint_line(ptx1, pty1, dwg)

				ptx2 = self.A10
				pty2 = ptx2 + (100) * np.array([0, -1])
				self.data_object.draw_faint_line(ptx2, pty2, dwg)

			else:

				ptx1 = self.A1
				pty1 = ptx1 + (100) * np.array([0, -1])
				self.data_object.draw_faint_line(ptx1, pty1, dwg)

				ptx2 = self.A9
				pty2 = ptx2 + (100) * np.array([0, -1])
				self.data_object.draw_faint_line(ptx2, pty2, dwg)

			if self.data_object.conn_loc == "Flange" or self.data_object.conn_loc == "Back to Back Web":
				point1 = ptx2 - (self.data_object.weld_inline/4) * np.array([1, 0])
				params = {"offset": -100, "textoffset": 10, "lineori": "right", "endlinedim": 10, "arrowlen": 20}
				self.data_object.draw_dimension_outer_arrow(dwg,point1,ptx2, str(self.data_object.weld_inline/4),params)
			else:
				point1 = ptx2 - (self.data_object.weld_inline / 2) * np.array([1, 0])
				params = {"offset": -100, "textoffset": 10, "lineori": "right", "endlinedim": 10, "arrowlen": 20}
				self.data_object.draw_dimension_outer_arrow(dwg, point1, ptx2, str(self.data_object.weld_inline/2),
															params)

		# ------------------------------------------  Sectional arrow -------------------------------------------
		pt_a1 = self.A4 - (200) * np.array([0, -1]) - (100 * np.array([1, 0]))
		pt_b1 = pt_a1 + (50 * np.array([0, -1]))
		txt_1 = pt_b1 + (50 * np.array([0, -1])) + (20 * np.array([-1, 0]))
		text = "C"
		self.data_object.draw_cross_section(dwg, pt_a1, pt_b1, txt_1, text)

		pt_a2 = self.A3 - (200) * np.array([0, -1]) + (100 * np.array([1, 0]))
		pt_b2 = pt_a2 + (50 * np.array([0, -1]))
		txt_2 = pt_b2 + (50 * np.array([0, -1])) + (20 * np.array([-1, 0]))
		self.data_object.draw_cross_section(dwg, pt_a2, pt_b2, txt_2, text)

		dwg.add(dwg.line(pt_a1, pt_a2).stroke('black', width=1.5, linecap='square'))

		if self.data_object.conn_loc == "Back to Back Web":

			pt_a3 = self.A3 + (100) * np.array([0, 1]) + (300 * np.array([1, 0]))
			pt_b3 = pt_a3 + (50 * np.array([-1, 0]))
			txt_3 = pt_b3 + (75 * np.array([-1, 0])) + (20 * np.array([0, 1]))
			text = "B"
			self.data_object.draw_cross_section(dwg, pt_a3, pt_b3, txt_3, text)

			pt_a4 = self.A16 - (100) * np.array([0, 1]) + (300 * np.array([1, 0]))
			pt_b4 = pt_a4 + (50 * np.array([-1, 0]))
			txt_4 = pt_b4 + (75 * np.array([-1, 0])) + (20 * np.array([0, 1]))
			self.data_object.draw_cross_section(dwg, pt_a4, pt_b4, txt_4, text)

			dwg.add(dwg.line(pt_a3, pt_a4).stroke('black', width=1.5, linecap='square'))
		else:
			pt_a3 = self.A2 - (100) * np.array([0, 1]) + (300 * np.array([1, 0]))
			pt_b3 = pt_a3 + (50 * np.array([-1, 0]))
			txt_3 = pt_b3 + (75 * np.array([-1, 0])) + (20 * np.array([0, 1]))
			text = "B"
			self.data_object.draw_cross_section(dwg, pt_a3, pt_b3, txt_3, text)

			pt_a4 = self.A3 + (100) * np.array([0, 1]) + (300 * np.array([1, 0]))
			pt_b4 = pt_a4 + (50 * np.array([-1, 0]))
			txt_4 = pt_b4 + (75 * np.array([-1, 0])) + (20 * np.array([0, 1]))
			self.data_object.draw_cross_section(dwg, pt_a4, pt_b4, txt_4, text)

			dwg.add(dwg.line(pt_a3, pt_a4).stroke('black', width=1.5, linecap='square'))

		# ------------------------------------------  View details -------------------------------------------
		ptx = self.A4 + (self.data_object.member_length/4)* np.array([1, 0]) + 300 * np.array([0, 1])
		dwg.add(dwg.text('Top view (Sec A-A) ', insert=ptx, fill='black', font_family="sans-serif", font_size=28))
		ptx1 = ptx + 40 * np.array([0, 1])
		dwg.add(dwg.text('(All dimensions are in "mm")', insert=ptx1, fill='black', font_family="sans-serif", font_size=28))
#
		dwg.save()


#
#
class Side_View (object):
	"""
	Contains functions for generating the side view of the Extended bothway endplate connection.

	"""

	def __init__(self, extnd_common_object):
		self.data_object = extnd_common_object

		# =========================  End Plate 1  =========================
		if self.data_object.conn_loc == "Leg":
			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			if self.data_object.weld_oppline > self.data_object.leg_min:
				ptA2x = ptA1x
				ptA2y = ptA1y + self.data_object.leg_max
				self.A2 = np.array([ptA2x, ptA2y])

				ptA3x = ptA2x + self.data_object.leg_min
				ptA3y = ptA2y
				self.A3 = np.array([ptA3x, ptA3y])

				ptA4x = ptA3x
				ptA4y = ptA3y - self.data_object.t
				self.A4 = np.array([ptA4x, ptA4y])

				ptA5x = ptA4x - self.data_object.leg_min + self.data_object.t
				ptA5y = ptA4y
				self.A5 = np.array([ptA5x, ptA5y])

				ptA6x = ptA1x + self.data_object.t
				ptA6y = ptA1y
				self.A6 = np.array([ptA6x, ptA6y])

				ptA7x = ptA1x
				ptA7y = ptA1y + self.data_object.leg_max/2 - 2 * self.data_object.weld_oppline
				self.A7 = np.array([ptA7x, ptA7y])

				ptA8x = ptA2x
				ptA8y = ptA2y - self.data_object.leg_max/2 + 2 * self.data_object.weld_oppline
				self.A8 = np.array([ptA8x, ptA8y])

				ptA9x = ptA8x - self.data_object.t
				ptA9y = ptA8y
				self.A9 = np.array([ptA9x, ptA9y])

				ptA10x = ptA7x - self.data_object.t
				ptA10y = ptA7y
				self.A10 = np.array([ptA10x, ptA10y])


			else:
				ptA2x = ptA1x
				ptA2y = ptA1y + self.data_object.leg_min
				self.A2 = np.array([ptA2x, ptA2y])

				ptA3x = ptA2x + self.data_object.leg_max
				ptA3y = ptA2y
				self.A3 = np.array([ptA3x, ptA3y])

				ptA4x = ptA3x
				ptA4y = ptA3y - self.data_object.t
				self.A4 = np.array([ptA4x, ptA4y])

				ptA5x = ptA4x - self.data_object.leg_max + self.data_object.t
				ptA5y = ptA4y
				self.A5 = np.array([ptA5x, ptA5y])

				ptA6x = ptA5x
				ptA6y = ptA1y
				self.A6 = np.array([ptA6x, ptA6y])

				ptA7x = ptA1x
				ptA7y = ptA1y + self.data_object.leg_min/2 - 2 * self.data_object.weld_oppline
				self.A7 = np.array([ptA7x, ptA7y])

				ptA8x = ptA2x
				ptA8y = ptA2y - self.data_object.leg_min/2 + 2 * self.data_object.weld_oppline
				self.A8 = np.array([ptA8x, ptA8y])

				ptA9x = ptA8x - self.data_object.t
				ptA9y = ptA8y
				self.A9 = np.array([ptA9x, ptA9y])

				ptA10x = ptA7x - self.data_object.t
				ptA10y = ptA7y
				self.A10 = np.array([ptA10x, ptA10y])

		elif self.data_object.conn_loc == "Back to Back Angles":
			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			if self.data_object.weld_oppline > self.data_object.leg_min:
				ptA2x = ptA1x
				ptA2y = ptA1y + self.data_object.leg_max
				self.A2 = np.array([ptA2x, ptA2y])

				ptA3x = ptA2x + self.data_object.leg_min
				ptA3y = ptA2y
				self.A3 = np.array([ptA3x, ptA3y])

				ptA4x = ptA3x
				ptA4y = ptA3y - self.data_object.t
				self.A4 = np.array([ptA4x, ptA4y])

				ptA5x = ptA4x - self.data_object.leg_min + self.data_object.t
				ptA5y = ptA4y
				self.A5 = np.array([ptA5x, ptA5y])

				ptA6x = ptA5x
				ptA6y = ptA1y
				self.A6 = np.array([ptA6x, ptA6y])

				ptA7x = ptA1x
				ptA7y = ptA1y + self.data_object.leg_max/2 - self.data_object.weld_oppline
				self.A7 = np.array([ptA7x, ptA7y])

				ptA8x = ptA2x
				ptA8y = ptA2y - self.data_object.leg_max/2 + self.data_object.weld_oppline
				self.A8 = np.array([ptA8x, ptA8y])

				ptA9x = ptA8x - self.data_object.t
				ptA9y = ptA8y
				self.A9 = np.array([ptA9x, ptA9y])

				ptA10x = ptA7x - self.data_object.t
				ptA10y = ptA7y
				self.A10 = np.array([ptA10x, ptA10y])

				ptA11x = ptA1x - self.data_object.t
				ptA11y = ptA1y
				self.A11 = np.array([ptA11x, ptA11y])

				ptA12x = ptA2x - self.data_object.t
				ptA12y = ptA2y
				self.A12 = np.array([ptA12x, ptA12y])

				ptA13x = ptA12x - self.data_object.leg_min
				ptA13y = ptA12y
				self.A13 = np.array([ptA13x, ptA13y])

				ptA14x = ptA13x
				ptA14y = ptA13y - self.data_object.t
				self.A14 = np.array([ptA14x, ptA14y])

				ptA15x = ptA14x + self.data_object.leg_min - self.data_object.t
				ptA15y = ptA14y
				self.A15 = np.array([ptA15x, ptA15y])

				ptA16x = ptA15x
				ptA16y = ptA11y
				self.A16 = np.array([ptA16x, ptA16y])

			else:
				ptA2x = ptA1x
				ptA2y = ptA1y + self.data_object.leg_min
				self.A2 = np.array([ptA2x, ptA2y])

				ptA3x = ptA2x + self.data_object.leg_max
				ptA3y = ptA2y
				self.A3 = np.array([ptA3x, ptA3y])

				ptA4x = ptA3x
				ptA4y = ptA3y - self.data_object.t
				self.A4 = np.array([ptA4x, ptA4y])

				ptA5x = ptA4x - self.data_object.leg_max + self.data_object.t
				ptA5y = ptA4y
				self.A5 = np.array([ptA5x, ptA5y])

				ptA6x = ptA5x
				ptA6y = ptA1y
				self.A6 = np.array([ptA6x, ptA6y])

				ptA7x = ptA1x
				ptA7y = ptA1y + self.data_object.leg_min/2 -  self.data_object.weld_oppline
				self.A7 = np.array([ptA7x, ptA7y])

				ptA8x = ptA2x
				ptA8y = ptA2y - self.data_object.leg_min/2 + self.data_object.weld_oppline
				self.A8 = np.array([ptA8x, ptA8y])

				ptA9x = ptA8x - self.data_object.t
				ptA9y = ptA8y
				self.A9 = np.array([ptA9x, ptA9y])

				ptA10x = ptA7x - self.data_object.t
				ptA10y = ptA7y
				self.A10 = np.array([ptA10x, ptA10y])

				ptA11x = ptA1x - self.data_object.t
				ptA11y = ptA1y
				self.A11 = np.array([ptA11x, ptA11y])

				ptA12x = ptA2x - self.data_object.t
				ptA12y = ptA2y
				self.A12 = np.array([ptA12x, ptA12y])

				ptA13x = ptA12x - self.data_object.leg_max
				ptA13y = ptA12y
				self.A13 = np.array([ptA13x, ptA13y])

				ptA14x = ptA13x
				ptA14y = ptA13y - self.data_object.t
				self.A14 = np.array([ptA14x, ptA14y])

				ptA15x = ptA14x + self.data_object.leg_max - self.data_object.t
				ptA15y = ptA14y
				self.A15 = np.array([ptA15x, ptA15y])

				ptA16x = ptA15x
				ptA16y = ptA11y
				self.A16 = np.array([ptA16x, ptA16y])

		elif self.data_object.conn_loc == "Star Angles":
			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			if self.data_object.weld_oppline > 2 * self.data_object.leg_min:
				ptA2x = ptA1x
				ptA2y = ptA1y + self.data_object.leg_max
				self.A2 = np.array([ptA2x, ptA2y])

				ptA3x = ptA2x + self.data_object.leg_min
				ptA3y = ptA2y
				self.A3 = np.array([ptA3x, ptA3y])

				ptA4x = ptA3x
				ptA4y = ptA3y - self.data_object.t
				self.A4 = np.array([ptA4x, ptA4y])

				ptA5x = ptA4x - self.data_object.leg_min + self.data_object.t
				ptA5y = ptA4y
				self.A5 = np.array([ptA5x, ptA5y])

				ptA6x = ptA5x
				ptA6y = ptA1y
				self.A6 = np.array([ptA6x, ptA6y])

				ptA7x = ptA1x
				ptA7y = ptA1y + self.data_object.leg_max/2 - self.data_object.weld_oppline
				self.A7 = np.array([ptA7x, ptA7y])

				ptA8x = ptA7x
				ptA8y = ptA7y + 2 * self.data_object.weld_oppline + self.data_object.leg_max
				self.A8 = np.array([ptA8x, ptA8y])

				ptA9x = ptA8x - self.data_object.t
				ptA9y = ptA8y
				self.A9 = np.array([ptA9x, ptA9y])

				ptA10x = ptA7x - self.data_object.t
				ptA10y = ptA7y
				self.A10 = np.array([ptA10x, ptA10y])

				ptA11x = ptA2x - self.data_object.t
				ptA11y = ptA2y
				self.A11 = np.array([ptA11x, ptA11y])

				ptA12x = ptA11x
				ptA12y = ptA11y + self.data_object.leg_max
				self.A12 = np.array([ptA12x, ptA12y])

				ptA13x = ptA12x - self.data_object.t
				ptA13y = ptA12y
				self.A13 = np.array([ptA13x, ptA13y])

				ptA14x = ptA13x
				ptA14y = ptA13y + self.data_object.t - self.data_object.leg_max
				self.A14 = np.array([ptA14x, ptA14y])

				ptA15x = ptA14x - self.data_object.leg_min + self.data_object.t
				ptA15y = ptA14y
				self.A15 = np.array([ptA15x, ptA15y])

				ptA16x = ptA15x
				ptA16y = ptA11y
				self.A16 = np.array([ptA16x, ptA16y])

			else:
				ptA2x = ptA1x
				ptA2y = ptA1y + self.data_object.leg_min
				self.A2 = np.array([ptA2x, ptA2y])

				ptA3x = ptA2x + self.data_object.leg_max
				ptA3y = ptA2y
				self.A3 = np.array([ptA3x, ptA3y])

				ptA4x = ptA3x
				ptA4y = ptA3y - self.data_object.t
				self.A4 = np.array([ptA4x, ptA4y])

				ptA5x = ptA4x - self.data_object.leg_max + self.data_object.t
				ptA5y = ptA4y
				self.A5 = np.array([ptA5x, ptA5y])

				ptA6x = ptA5x
				ptA6y = ptA1y
				self.A6 = np.array([ptA6x, ptA6y])

				ptA7x = ptA1x
				ptA7y = ptA1y + self.data_object.leg_min/2 - self.data_object.weld_oppline
				self.A7 = np.array([ptA7x, ptA7y])

				ptA8x = ptA7x
				ptA8y = ptA7y + 2 * self.data_object.weld_oppline + self.data_object.leg_min
				self.A8 = np.array([ptA8x, ptA8y])

				ptA9x = ptA8x - self.data_object.t
				ptA9y = ptA8y
				self.A9 = np.array([ptA9x, ptA9y])

				ptA10x = ptA7x - self.data_object.t
				ptA10y = ptA7y
				self.A10 = np.array([ptA10x, ptA10y])

				ptA11x = ptA2x - self.data_object.t
				ptA11y = ptA2y
				self.A11 = np.array([ptA11x, ptA11y])

				ptA12x = ptA11x
				ptA12y = ptA11y + self.data_object.leg_min
				self.A12 = np.array([ptA12x, ptA12y])

				ptA13x = ptA12x - self.data_object.t
				ptA13y = ptA12y
				self.A13 = np.array([ptA13x, ptA13y])

				ptA14x = ptA13x
				ptA14y = ptA13y + self.data_object.t - self.data_object.leg_min
				self.A14 = np.array([ptA14x, ptA14y])

				ptA15x = ptA14x - self.data_object.leg_max + self.data_object.t
				ptA15y = ptA14y
				self.A15 = np.array([ptA15x, ptA15y])

				ptA16x = ptA15x
				ptA16y = ptA11y
				self.A16 = np.array([ptA16x, ptA16y])

		elif self.data_object.section_type == "Channels" and (self.data_object.conn_loc == "Back to Back Web" or self.data_object.conn_loc == "Web"):
			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			ptA2x = ptA1x
			ptA2y = ptA1y + self.data_object.member_d
			self.A2 = np.array([ptA2x, ptA2y])

			ptA3x = ptA2x + self.data_object.member_B
			ptA3y = ptA2y
			self.A3 = np.array([ptA3x, ptA3y])

			ptA4x = ptA3x
			ptA4y = ptA3y - self.data_object.member_tf
			self.A4 = np.array([ptA4x, ptA4y])

			ptA5x = ptA4x - self.data_object.member_B + self.data_object.member_tw
			ptA5y = ptA4y
			self.A5 = np.array([ptA5x, ptA5y])

			ptA6x = ptA5x
			ptA6y = ptA5y - self.data_object.member_d + 2 * self.data_object.member_tf
			self.A6 = np.array([ptA6x, ptA6y])

			ptA7x = ptA6x + self.data_object.member_B - self.data_object.member_tw
			ptA7y = ptA6y
			self.A7 = np.array([ptA7x, ptA7y])

			ptA8x = ptA7x
			ptA8y = ptA7y - self.data_object.member_tf
			self.A8 = np.array([ptA8x, ptA8y])

			ptA9x = ptA1x
			ptA9y = ptA1y - self.data_object.member_d/2
			self.A9 = np.array([ptA9x, ptA9y])

			ptA10x = ptA2x
			ptA10y = ptA2y + self.data_object.member_d/2
			self.A10 = np.array([ptA10x, ptA10y])

			if self.data_object.section_type == "Channels":

				ptA11x = ptA10x - self.data_object.member_tw
				ptA11y = ptA10y
				self.A11 = np.array([ptA11x, ptA11y])

				ptA12x = ptA9x - self.data_object.member_tw
				ptA12y = ptA9y
				self.A12 = np.array([ptA12x, ptA12y])
			else:
				ptA11x = ptA10x - self.data_object.plate_thickness
				ptA11y = ptA10y
				self.A11 = np.array([ptA11x, ptA11y])

				ptA12x = ptA9x - self.data_object.plate_thickness
				ptA12y = ptA9y
				self.A12 = np.array([ptA12x, ptA12y])

			if  self.data_object.conn_loc == "Back to Back Web":
				ptA13x = ptA1x - self.data_object.plate_thickness
				ptA13y = ptA1y
				self.A13 = np.array([ptA13x, ptA13y])

				ptA14x = ptA2x - self.data_object.plate_thickness
				ptA14y = ptA2y
				self.A14 = np.array([ptA14x, ptA14y])

				ptA15x = ptA14x - self.data_object.member_B
				ptA15y = ptA14y
				self.A15 = np.array([ptA15x, ptA15y])

				ptA16x = ptA15x
				ptA16y = ptA15y - self.data_object.member_tf
				self.A16 = np.array([ptA16x, ptA16y])

				ptA17x = ptA16x + self.data_object.member_B - self.data_object.member_tw
				ptA17y = ptA16y
				self.A17 = np.array([ptA17x, ptA17y])

				ptA18x = ptA17x
				ptA18y = ptA17y - self.data_object.member_d + 2* self.data_object.member_tf
				self.A18 = np.array([ptA18x, ptA18y])

				ptA19x = ptA18x - self.data_object.member_B + self.data_object.member_tw
				ptA19y = ptA18y
				self.A19 = np.array([ptA19x, ptA19y])

				ptA20x = ptA19x
				ptA20y = ptA19y - self.data_object.member_tf
				self.A20 = np.array([ptA20x, ptA20y])


		elif self.data_object.section_type == "Channels" and self.data_object.conn_loc =="Flange":

			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			ptA2x = ptA1x
			ptA2y = ptA1y + self.data_object.member_d
			self.A2 = np.array([ptA2x, ptA2y])

			ptA3x = ptA2x + self.data_object.member_B
			ptA3y = ptA2y
			self.A3 = np.array([ptA3x, ptA3y])

			ptA4x = ptA3x
			ptA4y = ptA3y - self.data_object.member_tf
			self.A4 = np.array([ptA4x, ptA4y])

			ptA5x = ptA4x - self.data_object.member_B + self.data_object.member_tw
			ptA5y = ptA4y
			self.A5 = np.array([ptA5x, ptA5y])

			ptA6x = ptA5x
			ptA6y = ptA5y - self.data_object.member_d + 2 * self.data_object.member_tf
			self.A6 = np.array([ptA6x, ptA6y])

			ptA7x = ptA6x + self.data_object.member_B - self.data_object.member_tw
			ptA7y = ptA6y
			self.A7 = np.array([ptA7x, ptA7y])

			ptA8x = ptA7x
			ptA8y = ptA7y - self.data_object.member_tf
			self.A8 = np.array([ptA8x, ptA8y])

			ptA9x = ptA2x + self.data_object.member_B/2 - self.data_object.weld_oppline/4
			ptA9y = ptA2y
			self.A9 = np.array([ptA9x, ptA9y])

			ptA10x = ptA9x
			ptA10y = ptA9y + self.data_object.member_tf
			self.A10 = np.array([ptA10x, ptA10y])

			ptA11x = ptA10x + self.data_object.weld_oppline/2
			ptA11y = ptA10y
			self.A11 = np.array([ptA11x, ptA11y])

			ptA12x = ptA11x
			ptA12y = ptA9y
			self.A12 = np.array([ptA12x, ptA12y])

			ptA13x = ptA1x + self.data_object.member_B / 2 - self.data_object.weld_oppline /4
			ptA13y = ptA1y
			self.A13 = np.array([ptA13x, ptA13y])

			ptA14x = ptA13x
			ptA14y = ptA13y - self.data_object.member_tf
			self.A14 = np.array([ptA14x, ptA14y])

			ptA15x = ptA14x + self.data_object.weld_oppline / 2
			ptA15y = ptA14y
			self.A15 = np.array([ptA15x, ptA15y])

			ptA16x = ptA15x
			ptA16y = ptA13y
			self.A16 = np.array([ptA16x, ptA16y])

		elif self.data_object.section_type == "Beams" or "Columns":

			ptA1x = 0
			ptA1y = 0
			self.A1 = np.array([ptA1x, ptA1y])

			ptA2x = ptA1x + self.data_object.member_B
			ptA2y = ptA1y
			self.A2 = np.array([ptA2x, ptA2y])

			ptA3x = ptA2x
			ptA3y = ptA2y + self.data_object.member_tf
			self.A3 = np.array([ptA3x, ptA3y])

			ptA4x = ptA1x
			ptA4y = ptA3y
			self.A4 = np.array([ptA4x, ptA4y])

			ptA5x = ptA4x + self.data_object.member_B/2 - self.data_object.member_tw/2
			ptA5y = ptA4y
			self.A5 = np.array([ptA5x, ptA5y])

			ptA6x = ptA5x + self.data_object.member_tw
			ptA6y = ptA5y
			self.A6 = np.array([ptA6x, ptA6y])

			ptA7x = ptA6x
			ptA7y = ptA6y + self.data_object.member_d - 2 * self.data_object.member_tf
			self.A7 = np.array([ptA7x, ptA7y])

			ptA8x = ptA7x - self.data_object.member_tw
			ptA8y = ptA7y
			self.A8 = np.array([ptA8x, ptA8y])

			ptA9x = ptA1x
			ptA9y = ptA1y + self.data_object.member_d - self.data_object.member_tf
			self.A9 = np.array([ptA9x, ptA9y])

			ptA10x = ptA9x + self.data_object.member_B
			ptA10y = ptA9y
			self.A10 = np.array([ptA10x, ptA10y])

			ptA11x = ptA10x
			ptA11y = ptA10y + self.data_object.member_tf
			self.A11 = np.array([ptA11x, ptA11y])

			ptA12x = ptA9x
			ptA12y = ptA11y
			self.A12 = np.array([ptA12x, ptA12y])

			if self.data_object.conn_loc =="Flange":

				ptA13x = ptA1x + self.data_object.member_B/2 - self.data_object.weld_oppline/4
				ptA13y = ptA1y
				self.A13 = np.array([ptA13x, ptA13y])

				ptA14x = ptA13x
				ptA14y = ptA13y - self.data_object.member_tf
				self.A14 = np.array([ptA14x, ptA14y])

				ptA15x = ptA14x + self.data_object.weld_oppline / 2
				ptA15y = ptA14y
				self.A15 = np.array([ptA15x, ptA15y])

				ptA16x = ptA15x
				ptA16y = ptA13y
				self.A16 = np.array([ptA16x, ptA16y])

				ptA17x = ptA12x + self.data_object.member_B / 2 - self.data_object.weld_oppline /4
				ptA17y = ptA12y
				self.A17 = np.array([ptA17x, ptA17y])

				ptA18x = ptA17x
				ptA18y = ptA17y + self.data_object.member_tf
				self.A18 = np.array([ptA18x, ptA18y])

				ptA19x = ptA18x + self.data_object.weld_oppline / 2
				ptA19y = ptA18y
				self.A19 = np.array([ptA19x, ptA19y])

				ptA20x = ptA19x
				ptA20y = ptA17y
				self.A20 = np.array([ptA20x, ptA20y])

			else:
				ptA13x = ptA6x
				ptA13y = ptA6y + self.data_object.member_d/2 - self.data_object.member_tf - self.data_object.weld_oppline/2
				self.A13 = np.array([ptA13x, ptA13y])

				ptA14x = ptA13x + self.data_object.member_tw
				ptA14y = ptA13y
				self.A14 = np.array([ptA14x, ptA14y])

				ptA15x = ptA14x
				ptA15y = ptA14y + self.data_object.weld_oppline
				self.A15 = np.array([ptA15x, ptA15y])

				ptA16x = ptA13x
				ptA16y = ptA15y
				self.A16 = np.array([ptA16x, ptA16y])

		else:
			pass

#
	def call_Side_View(self, filename):
		"""

		Args:
			filename: path of the images to be saved

		Returns:
			Saves the image in the folder

		"""
		if self.data_object.conn_loc == "Leg":
			# wd = int((self.data_object.member_length) + 1000)
			if self.data_object.weld_oppline > self.data_object.leg_min:
				ht = int(self.data_object.leg_max + 400)
				wd = int(self.data_object.leg_max + 700)
			else:
				ht = int(self.data_object.leg_min + 400)
				wd = int(self.data_object.leg_min + 700)
			dwg = svgwrite.Drawing(filename, size=('100%', '100%'), viewBox=(
				'-350 -200 {} {}').format(wd, ht))
			dwg.add(dwg.polyline(points=[self.A1, self.A2, self.A3, self.A4, self.A5, self.A6,self.A1],
				stroke='blue', fill="none", stroke_width=1))
			# dwg.add(dwg.line(self.A6, self.A7).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A7, self.A1).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A2, self.A8).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A8, self.A9).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A9, self.A10).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A7, self.A10).stroke('blue', width=1, linecap='square'))
			# dwg.add(dwg.line(self.A1, self.A4).stroke('blue', width=1, linecap='square'))
		elif self.data_object.conn_loc == "Back to Back Angles":
			# wd = int((self.data_object.member_length) + 1000)
			if self.data_object.weld_oppline > self.data_object.leg_min:
				ht = int(self.data_object.leg_max + 500)
				wd = int(self.data_object.leg_max + 750)
			else:
				ht = int(self.data_object.leg_min + 500)
				wd = int(self.data_object.leg_min + 750)
			dwg = svgwrite.Drawing(filename, size=('100%', '100%'), viewBox=(
				'-375 -250 {} {}').format(wd, ht))
			dwg.add(dwg.polyline(points=[self.A1, self.A2, self.A3, self.A4, self.A5, self.A6,self.A1],
				stroke='blue', fill="none", stroke_width=1))
			dwg.add(dwg.polyline(points=[self.A11, self.A12, self.A13, self.A14, self.A15, self.A16, self.A11],
								 stroke='blue', fill="none", stroke_width=1))
			# dwg.add(dwg.line(self.A6, self.A7).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A7, self.A1).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A2, self.A8).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A9, self.A8).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A12, self.A9).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A11, self.A10).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A7, self.A10).stroke('blue', width=1, linecap='square'))

		elif self.data_object.conn_loc == "Star Angles":
			# wd = int((self.data_object.membe600r_length) + 1000)
			if self.data_object.weld_oppline > self.data_object.leg_min:
				ht = int(2*self.data_object.leg_max + 500)
				wd = int(self.data_object.leg_max + 750)
			else:
				ht = int(2*self.data_object.leg_min + 500)
				wd = int(self.data_object.leg_min + 750)
			dwg = svgwrite.Drawing(filename, size=('100%', '100%'), viewBox=(
				'-375 -250 {} {}').format(wd, ht))
			dwg.add(dwg.polyline(points=[self.A1, self.A2, self.A3, self.A4, self.A5, self.A6,self.A1],
				stroke='blue', fill="none", stroke_width=1))
			dwg.add(dwg.polyline(points=[self.A11, self.A12, self.A13, self.A14, self.A15, self.A16, self.A11],
								 stroke='blue', fill="none", stroke_width=1))
			# dwg.add(dwg.line(self.A6, self.A7).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A7, self.A1).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A2, self.A8).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A9, self.A8).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A12, self.A9).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A11, self.A10).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A7, self.A10).stroke('blue', width=1, linecap='square'))

		elif self.data_object.section_type == "Channels":
			# wd = int((self.data_object.member_length) + 1000)

			ht = int(self.data_object.member_d + 600)
			wd = int(self.data_object.member_B + 1000)
			dwg = svgwrite.Drawing(filename, size=('100%', '100%'), viewBox=(
				'-500 -300 {} {}').format(wd, ht))
			dwg.add(dwg.polyline(points=[self.A1, self.A2, self.A3, self.A4, self.A5, self.A6,self.A7,self.A8,self.A1],
				stroke='blue', fill="none", stroke_width=1))
			if self.data_object.conn_loc == "Back to Back Web" or self.data_object.conn_loc == "Web":
				if self.data_object.conn_loc == "Back to Back Web":
					dwg.add(dwg.line(self.A13, self.A20).stroke('blue', width=1, linecap='square'))
					dwg.add(dwg.line(self.A19, self.A20).stroke('blue', width=1, linecap='square'))
					dwg.add(dwg.line(self.A19, self.A18).stroke('blue', width=1, linecap='square'))
					dwg.add(dwg.line(self.A17, self.A18).stroke('blue', width=1, linecap='square'))
					dwg.add(dwg.line(self.A17, self.A16).stroke('blue', width=1, linecap='square'))
					dwg.add(dwg.line(self.A15, self.A16).stroke('blue', width=1, linecap='square'))
					dwg.add(dwg.line(self.A15, self.A14).stroke('blue', width=1, linecap='square'))

				dwg.add(dwg.line(self.A9, self.A1).stroke('blue', width=1, linecap='square'))
				dwg.add(dwg.line(self.A2, self.A10).stroke('blue', width=1, linecap='square'))
				dwg.add(dwg.line(self.A10, self.A11).stroke('blue', width=1, linecap='square'))
				dwg.add(dwg.line(self.A12, self.A11).stroke('blue', width=1, linecap='square'))
				dwg.add(dwg.line(self.A12, self.A9).stroke('blue', width=1, linecap='square'))
			else:
				dwg.add(dwg.line(self.A9, self.A10).stroke('blue', width=1, linecap='square'))
				dwg.add(dwg.line(self.A11, self.A10).stroke('blue', width=1, linecap='square'))
				dwg.add(dwg.line(self.A12, self.A11).stroke('blue', width=1, linecap='square'))
				dwg.add(dwg.line(self.A13, self.A14).stroke('blue', width=1, linecap='square'))
				dwg.add(dwg.line(self.A14, self.A15).stroke('blue', width=1, linecap='square'))
				dwg.add(dwg.line(self.A15, self.A16).stroke('blue', width=1, linecap='square'))

		elif self.data_object.section_type == "Beams" or "Columns":
			ht = int(self.data_object.member_d + 600)
			wd = int(self.data_object.member_B + 1000)
			dwg = svgwrite.Drawing(filename, size=('100%', '100%'), viewBox= ('-500 -300 {} {}').format(wd, ht))
			dwg.add(dwg.polyline(points=[self.A1, self.A2, self.A3, self.A4, self.A1], stroke='blue', fill="none", stroke_width=1))
			dwg.add(dwg.polyline(points=[self.A9, self.A10, self.A11, self.A12, self.A9], stroke='blue', fill="none", stroke_width=1))
			dwg.add(dwg.line(self.A5, self.A8).stroke('blue', width=1, linecap='square'))
			dwg.add(dwg.line(self.A7, self.A6).stroke('blue', width=1, linecap='square'))
			if self.data_object.conn_loc == "Flange":
				dwg.add(dwg.line(self.A13, self.A14).stroke('blue', width=1, linecap='square'))
				dwg.add(dwg.line(self.A14, self.A15).stroke('blue', width=1, linecap='square'))
				dwg.add(dwg.line(self.A15, self.A16).stroke('blue', width=1, linecap='square'))
				dwg.add(dwg.line(self.A17, self.A18).stroke('blue', width=1, linecap='square'))
				dwg.add(dwg.line(self.A19, self.A18).stroke('blue', width=1, linecap='square'))
				dwg.add(dwg.line(self.A19, self.A20).stroke('blue', width=1, linecap='square'))
			else:
				dwg.add(dwg.line(self.A13, self.A14).stroke('blue', width=1, linecap='square'))
				dwg.add(dwg.line(self.A14, self.A15).stroke('blue', width=1, linecap='square'))
				dwg.add(dwg.line(self.A15, self.A16).stroke('blue', width=1, linecap='square'))
		else:
			pass

		# ------------------------------------------  Primary Beam 1& 2 -------------------------------------------
		if self.data_object.section_type == "Angles":
			point = self.A1
			theta = 60
			offset = 100
			textdown = " "
			textup =  str(self.data_object.beam_designation)
			element = " "
			self.data_object.draw_oriented_arrow_side(dwg, point, theta, "NW", offset, textup, textdown, element)
		else:
			point = self.A1
			theta = 60
			offset = 100
			textdown = " "
			textup = str(self.data_object.beam_designation)
			element = " "
			self.data_object.draw_oriented_arrow(dwg, point, theta, "NW", offset, textup, textdown, element)


		# ------------------------------------------Dimension Labelling -------------------------------------------
		if self.data_object.conn_loc == "Leg":
			ptx1 = self.A1
			pty1 = ptx1 + (100) * np.array([-1, 0])
			self.data_object.draw_faint_line(ptx1, pty1, dwg)

			ptx2 = self.A2
			pty2 = ptx2 + (100) * np.array([-1, 0])
			self.data_object.draw_faint_line(ptx2, pty2, dwg)

			point1 = ptx2 - (self.data_object.weld_oppline) * np.array([0, 1])
			params = {"offset": 100, "textoffset": 40, "lineori": "left", "endlinedim": 5, "arrowlen": 20}
			self.data_object.draw_dimension_outer_arrow_side(dwg,ptx2, point1,str(self.data_object.weld_oppline),params)

		elif self.data_object.conn_loc == "Back to Back Angles":
			ptx1 = self.A1
			pty1 = ptx1 + (100) * np.array([-1, 0])
			self.data_object.draw_faint_line(ptx1, pty1, dwg)

			ptx2 = self.A2
			pty2 = ptx2 + (100) * np.array([-1, 0])
			self.data_object.draw_faint_line(ptx2, pty2, dwg)

			point1 = ptx2 - (self.data_object.weld_oppline/2) * np.array([0, 1])
			params = {"offset": 100, "textoffset": 40, "lineori": "left", "endlinedim": 5, "arrowlen": 20}
			self.data_object.draw_dimension_outer_arrow_side(dwg,ptx2, point1,str(self.data_object.weld_oppline/2),params)

		elif self.data_object.conn_loc == "Star Angles":
			ptx1 = self.A1
			pty1 = ptx1 + (100) * np.array([-1, 0])
			self.data_object.draw_faint_line(ptx1, pty1, dwg)

			ptx2 = self.A12
			pty2 = ptx2 + (100) * np.array([-1, 0])
			self.data_object.draw_faint_line(ptx2, pty2, dwg)


			point1 = ptx2 - (self.data_object.weld_oppline) * np.array([0, 1])
			params = {"offset": 100, "textoffset":50, "lineori": "left", "endlinedim": 5, "arrowlen": 20}
			self.data_object.draw_dimension_outer_arrow_side(dwg,ptx2, point1,str(self.data_object.weld_oppline),params)

		elif self.data_object.conn_loc == "Back to Back Web":
			ptx1 = self.A20
			pty1 = ptx1 + (100) * np.array([-1, 0])
			self.data_object.draw_faint_line(ptx1, pty1, dwg)

			ptx2 = self.A15
			pty2 = ptx2 + (100) * np.array([-1, 0])
			self.data_object.draw_faint_line(ptx2, pty2, dwg)

			point1 = ptx2 - (self.data_object.weld_oppline / 2) * np.array([0, 1])
			params = {"offset": 100, "textoffset": 100, "lineori": "left", "endlinedim": 5, "arrowlen": 20}
			self.data_object.draw_dimension_outer_arrow_side(dwg, ptx2, point1, str(self.data_object.weld_oppline / 2),
														params)

		elif self.data_object.conn_loc == "Web" and self.data_object.section_type == "Channels":
			ptx1 = self.A1
			pty1 = ptx1 + (100) * np.array([-1, 0])
			self.data_object.draw_faint_line(ptx1, pty1, dwg)

			ptx2 = self.A2
			pty2 = ptx2 + (100) * np.array([-1, 0])
			self.data_object.draw_faint_line(ptx2, pty2, dwg)

			point1 = ptx2 - (self.data_object.weld_oppline) * np.array([0, 1])
			params = {"offset": 100, "textoffset": 75, "lineori": "left", "endlinedim": 5, "arrowlen": 20}
			self.data_object.draw_dimension_outer_arrow_side(dwg, ptx2, point1, str(self.data_object.weld_oppline),
														params)

		elif self.data_object.conn_loc == "Flange":
			ptx1 = self.A14
			pty1 = ptx1 + (100) * np.array([0, -1])
			self.data_object.draw_faint_line(ptx1, pty1, dwg)

			ptx2 = self.A15
			pty2 = ptx2 + (100) * np.array([0, -1])
			self.data_object.draw_faint_line(ptx2, pty2, dwg)

			point1 = ptx2 - (self.data_object.weld_oppline/2) * np.array([1, 0])
			params = {"offset": 100, "textoffset": 40, "lineori": "left", "endlinedim": 5, "arrowlen": 20}
			self.data_object.draw_dimension_outer_arrow(dwg, point1, ptx2,str(self.data_object.weld_oppline/2),
														params)

		elif self.data_object.conn_loc == "Web" and (self.data_object.section_type == "Beams" or self.data_object.section_type =="Columns"):
			ptx1 = self.A13
			pty1 = ptx1 + (100) * np.array([-1, 0])
			self.data_object.draw_faint_line(ptx1, pty1, dwg)

			ptx2 = self.A16
			pty2 = ptx2 + (100) * np.array([-1, 0])
			self.data_object.draw_faint_line(ptx2, pty2, dwg)

			point1 = ptx2 - (self.data_object.weld_oppline) * np.array([0, 1])
			params = {"offset": 100, "textoffset": 100, "lineori": "left", "endlinedim": 5, "arrowlen": 20}
			self.data_object.draw_dimension_outer_arrow_side(dwg,ptx2,point1,str(self.data_object.weld_oppline),
														params)
		else:
			pass


## ------------------------------------------  View details--------------------------------------------------------
		if self.data_object.section_type == "Angles":
			ptx = self.A8 + 50 * np.array([0, 1]) - (len('(All dimensions are in "mm")')) * np.array([1, 0])
			dwg.add(dwg.text('Side view (Sec B-B) ', insert=ptx, fill='black', font_family="sans-serif", font_size=15))
			ptx1 = ptx + 40 * np.array([0, 1])
			dwg.add(dwg.text('(All dimensions are in "mm")', insert=ptx1, fill='black', font_family="sans-serif", font_size=15))
		elif self.data_object.section_type == "Channels":
			ptx = self.A10 + 100 * np.array([0, 1]) - (4  *len('(All dimensions are in "mm")')) * np.array([1, 0])
			dwg.add(dwg.text('Side view (Sec B-B) ', insert=ptx, fill='black', font_family="sans-serif", font_size=28))
			ptx1 = ptx + 40 * np.array([0, 1])
			dwg.add(dwg.text('(All dimensions are in "mm")', insert=ptx1, fill='black', font_family="sans-serif",
							 font_size=28))
		elif self.data_object.section_type == "Beams" or "Columns":
			ptx = self.A12 + 100 * np.array([0, 1]) - (len('(All dimensions are in "mm")')) * np.array([1, 0])
			dwg.add(dwg.text('Side view (Sec B-B) ', insert=ptx, fill='black', font_family="sans-serif", font_size=28))
			ptx1 = ptx + 40 * np.array([0, 1])
			dwg.add(dwg.text('(All dimensions are in "mm")', insert=ptx1, fill='black', font_family="sans-serif",
							 font_size=28))
		else:
			pass
#
		dwg.save()
