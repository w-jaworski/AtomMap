(*
 *  AtomMap: maps atoms in chemical reactions
 *  Copyright (C) 2015-2017 Wojciech Jaworski <wjaworski atSPAMfree mimuw dot edu dot pl> 
 *  Copyright (C) 2015-2017 Institute of Informatics, University of Warsaw                

 *  Copyright (C) 2015-2017 Sara Szymkuc <saraszymkuc atSPAMfree gmail dot com>          
 *  Copyright (C) 2015-2017 Barbara Mikulak <basia dot mikulak atSPAMfree gmail dot com> 
 *  Copyright (C) 2015-2017 Institute of Organic Chemistry Polish Academy of Sciences     

 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *)

open Types

let _ =
  if Array.length Sys.argv < 2 then print_endline "missing argument" else
  let query = Sys.argv.(1) in
  let re = {empty_record with rxn_id="id"; reaction_smile=query} in
  let messages,solutions = MatchingExec.map_atoms_pat re in
  let messages,solutions = if solutions = [] then MatchingExec.map_atoms false re else messages,solutions in
  let messages,solutions = if solutions = [] then MatchingExec.map_atoms true re else messages,solutions in
  let xml = Smiles.solutions_to_xml query messages solutions in
  print_endline (Xml.to_string_fmt xml)
