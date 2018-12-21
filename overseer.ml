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

let create_workers n_workers =
  let id = string_of_int (Unix.getpid ()) in
  let io_list = Int.fold 1 n_workers [] (fun io_list _ ->
    print_endline (id ^ " create_worker");
    let in_chan,out_chan = Unix.open_process "./worker" in
    let descr = Unix.descr_of_in_channel in_chan in
    (in_chan,out_chan,descr) :: io_list) in
  io_list

let io_list = create_workers 1

let execution io_list msg work (*file*) =
  let size = Xlist.size work in
  let r = ref (size (*+ Xlist.size io_list*)) in
  let size = string_of_int size in
  let work = ref work in
(*   let output = ref [] in *)
  let sum_result = ref [] in
  let id = string_of_int (Unix.getpid ()) in
  let descr_list = Xlist.map io_list (fun (_,_,descr) -> descr) in
  while !r <> 0 do
(*     print_endline (id ^ " Unix.select");  *)
    let list,_,_ = Unix.select descr_list [] [] (-1.) in
(*     print_endline (id ^ " selected " ^ (string_of_int (Xlist.size list)));  *)
    Xlist.iter list (fun descr2 ->
      decr r;
      Xlist.iter io_list (fun (in_chan,out_chan,descr) ->
        if descr = descr2 then (
          let idw = match Marshal.from_channel in_chan with
            Ready_to_work idw ->
(*               print_endline (idw ^ " ready"); *)
              idw
          | Work_done (idw,(messages,matchings)) ->
(*               print_endline (idw ^ " work done"); *)
              msg := messages @ !msg;
(*               output := s :: (!output); *)
(*               Exec.print_result file s; *)
              sum_result := matchings :: !sum_result;
(*               Exec.print_sum_result file !sum_result; *)
              idw in
          match !work with
            (id,params) :: l ->
              Marshal.to_channel out_chan (Work_with (id,params)) [Marshal.No_sharing];
              flush out_chan;
(*               print_endline (idw ^ " scheduled " ^ id ^ " of " ^ size); *)
              work := l
          | [] ->
(*               Marshal.to_channel out_chan Kill_yourself [Marshal.No_sharing];  *)
              (*print_endline (idw ^ " finished")*)())))
  done;
(*   print_endline (id ^ " exit"); *)
  !sum_result

(*let _ =
  if Array.length Sys.argv < 2 then print_endline "missing argument" else
  let n_workers = Paths.no_processes in
  let work = Exec.generate_queries Sys.argv.(1) Paths.lcg_timeout in
  File.file_out (Paths.results_path ^ "log") (fun file ->
  execution n_workers work file)*)
