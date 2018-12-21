(* Copyright (c) 2006-2017, Wojciech Jaworski
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 *)

open Types
open Printf

(*let _ =
(*   fprintf logfile "woo1\n%!"; *)
  let id = string_of_int (Unix.getpid ()) in
  let f = ref true in
(*   fprintf logfile "woo2 %s\n%!" id; *)
  Marshal.to_channel stdout (Ready_to_work id) [Marshal.No_sharing];
(*   fprintf logfile "woo3\n%!"; *)
  flush stdout;
(*   fprintf logfile "woo4\n%!"; *)
  while !f do
(*     fprintf logfile "woo5\n%!"; *)
    (match (try Marshal.from_channel stdin with e -> (fprintf logfile "%s\n%!" (Printexc.to_string e); Kill_yourself)) with
      Work_with (akt_id,(r,cand,r2,l)) ->
(*         print_endline (id ^ " working with " ^ akt_id); *)
(*         fprintf logfile "woo6\n%!"; *)
        let msg,matchings = (*try
          let matchings = CommonSubstructure.match_reaction 4 r2 labels in
          let matchings = AtomMapping.select_carbonyls r matchings in  (* !!! wyłączone testowo do OrgSyn *)
          [],Collection.flatten_map matchings (fun matching ->
            AtomMapping.match_hydrogens_and_fluors cand r2 matching)
          with e ->*) [(*Printexc.to_string e*)], Collection.empty in
(*         print_endline (id ^ " finished working with " ^ akt_id); *)
(*         fprintf logfile "woo7\n%!"; *)
        Marshal.to_channel stdout (Work_done(id, (msg,matchings))) [Marshal.No_sharing];
(*         fprintf logfile "woo8\n%!"; *)
        flush stdout
    | Kill_yourself -> f := false);
(*     fprintf logfile "woo9\n%!"; *)
  done;
(*   fprintf logfile "woo10\n%!"; *)
(*   print_endline (id ^ " work finished") *)
  ()*)

let _ =
(*   fprintf logfile "woo1\n%!"; *)
  let id = string_of_int (Unix.getpid ()) in
  let f = ref true in
(*   fprintf logfile "woo2 %s\n%!" id; *)
  Marshal.to_channel stdout (Ready_to_work id) [Marshal.No_sharing];
(*   fprintf logfile "woo3\n%!"; *)
  flush stdout;
(*   fprintf logfile "woo4\n%!"; *)
  while !f do
(*     fprintf logfile "woo5\n%!"; *)
    (match (try Marshal.from_channel stdin with e -> (fprintf logfile "%s\n%!" (Printexc.to_string e); Kill_yourself)) with
      Work_with (akt_id,(simple_flag,re)) ->
      Types.time := Sys.time ();
      Printf.fprintf logfile "A %s %s\n%!" re.rxn_id re.reaction_smile;
(*         print_endline (id ^ " working with " ^ akt_id); *)
(*         fprintf logfile "woo6\n%!"; *)
        let messages,solutions = try MatchingExec.map_atoms simple_flag re with e -> [Smiles.escape_string (Printexc.to_string e)],[] in
(*         print_endline (id ^ " finished working with " ^ akt_id); *)
(*         fprintf logfile "woo7\n%!"; *)
        Marshal.to_channel stdout (Work_done(id, (akt_id,re,messages,solutions))) [Marshal.No_sharing];
(*         fprintf logfile "woo8\n%!"; *)
      let time = Sys.time () -. !Types.time in
(*       Printf.fprintf logfile "B %s %s\nTIME: %f\n%!" re.rxn_id re.reaction_smile time; *)
      Printf.fprintf logfile "B %s %s\nTIME: %f\n%!" re.rxn_id re.reaction_smile time;
      Printf.fprintf logfile "STATS %s,%s,%f,%d,%s,%s\n%!" re.rxn_id re.reaction_smile time (Smiles.count_reactants_atoms re) (Import.string_of_is_correct re.is_correct) re.path;
        flush stdout
    | Kill_yourself -> f := false);
(*     fprintf logfile "woo9\n%!"; *)
  done;
(*   fprintf logfile "woo10\n%!"; *)
(*   print_endline (id ^ " work finished") *)
  ()
