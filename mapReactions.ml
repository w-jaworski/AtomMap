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

(*let _ =
  if Array.length Sys.argv < 4 then print_endline "missing argument\nusage: map_reactions name_prefix input_filename output_filename" else
  let reactions = Import.load_reactions Sys.argv.(2) in
  MatchingExec.map_atoms2 Sys.argv.(1) Sys.argv.(3) reactions*)

(*let _ =
  if Array.length Sys.argv < 5 then print_endline "missing argument\nusage: map_reactions no_processors name_prefix input_filename output_filename" else
  let reactions = Import.load_reactions Sys.argv.(3) in
  let no_processors = try int_of_string Sys.argv.(1) with _ -> failwith "invalid value for no_processors" in
  let simple_flag = false in
  MatchingExec.map_atoms2_distr no_processors simple_flag Sys.argv.(2) Sys.argv.(4) reactions*)

let _ =
  if Array.length Sys.argv < 3 then print_endline "missing argument\nusage: map_reactions no_processors input_filename" else
  let reactions = Import.load_reactions Sys.argv.(2) in
  let no_processors = try int_of_string Sys.argv.(1) with _ -> failwith "invalid value for no_processors" in
  let simple_flag = false in
  MatchingExec.map_atoms2_distr3 no_processors simple_flag reactions

(*let _ =
  if Array.length Sys.argv < 3 then print_endline "missing argument\nusage: map_reactions no_processors input_filename" else
  let reactions = Import.load_rxn2 Sys.argv.(2) in
  let no_processors = try int_of_string Sys.argv.(1) with _ -> failwith "invalid value for no_processors" in
  let simple_flag = false in
  MatchingExec.map_atoms2_distr3 no_processors simple_flag reactions*)

  
(*let _ =
  if Array.length Sys.argv < 3 then print_endline "missing argument\nusage: map_reactions no_processors input_filename" else
  let reactions = Import.load_wszystkie () in
  let no_processors = try int_of_string Sys.argv.(1) with _ -> failwith "invalid value for no_processors" in
  let simple_flag = false in
  MatchingExec.map_atoms2_distr3 no_processors simple_flag reactions*)
