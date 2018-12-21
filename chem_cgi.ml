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
open Xstd

let get_input () =
  let r = ref [] in
  (try
    while true do
      r := (input_line stdin) :: (!r)
    done;
    !r
  with End_of_file -> !r)

let rec translate_input_rec buf i size query =
  if i >= size then Buffer.contents buf else (
  let c,i =
    if String.get query i = '%' then
      Scanf.sscanf (String.sub query (i+1) 2) "%x" (fun a -> Char.chr a), i+3
    else String.get query i, i+1 in
  Buffer.add_char buf c;
  translate_input_rec buf i size query)

let translate_input query =
  match query with
    [query] ->
      if String.sub query 0 6 = "text0=" then
        let buf = Buffer.create (String.length query) in
        translate_input_rec buf 6 (String.length query) query
      else failwith "translate_input 1"
  | _ -> failwith "translate_input 2"

let generate_header () =
  Printf.printf "Content-type: text/xml\n";
  Printf.printf "Access-Control-Allow-Origin: *\n";
  Printf.printf "Access-Control-Request-Headers: x-requested-with\n";
  Printf.printf "Access-Control-Allow-Methods: POST, GET, OPTIONS\n";
  Printf.printf "\n"

let generate_trailer () =
  (*Printf.printf "</BODY>\n</HTML>\n"*)()

(* let host = "localhost" *)
(* let host = "wloczykij" *)
let host = "213.135.37.219"
let port = 2727

let get_sock_addr host_name port =
  let he = Unix.gethostbyname host_name in
  let addr = he.Unix.h_addr_list in
  Unix.ADDR_INET(addr.(0),port)

let _ =
  generate_header ();
  (try
    let query = get_input () in
    let query = translate_input query in
    let ic,oc = Unix.open_connection (get_sock_addr host port) in
    Printf.fprintf oc "%s\n%!" query;
    let messages,solutions = (Marshal.from_channel ic : (string list * (Types.reaction * Labels.t * int IntMap.t IntMap.t * int IntMap.t IntMap.t * int IntMap.t) list)) in
(*     Printf.fprintf oc "\n%!";  *)
    let _ = Unix.shutdown_connection ic in
    let xml = Smiles.solutions_to_xml query messages solutions in
    print_endline (Xml.to_string_fmt xml)
  with e -> print_endline (Printexc.to_string e));
  generate_trailer ()
