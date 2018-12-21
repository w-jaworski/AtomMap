(*
 *  AtomMap: maps atoms in chemical reactions
 *  Copyright (C) 2015-2017 Wojciech Jaworski <wjaworski atSPAMfree mimuw dot edu dot pl> 
 *  Copyright (C) 2015-2017 Institute of Informatics, University of Warsaw                

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


open Xstd
open Types

exception TreeToBig

(* Poniższe procedury szukają izomorfizmów zadanych podgrafów ignorując krotność wiązań i etykietując wierzchołki haszami *)

val are_isomorphic : bool -> int ->
           hash array ->
           hash array ->
           graph ->
           graph ->
           int Collection.collection -> int Collection.collection -> bool

val are_isomorphic_lists : bool -> int ->
          'a * (hash array * graph * int Collection.collection) list ->
          'a * (hash array * graph * int Collection.collection) list -> bool

val are_isomorphic_alt_lists : bool -> int -> 'a * (graph * hash array) -> 'a * (graph * hash array) ->
          (int Collection.collection * int Collection.collection) Collection.collection list ->
          bool

val find_isomorphisms : bool -> int ->
           hash array ->
           hash array ->
           graph ->
           graph ->
           int Collection.collection -> int Collection.collection -> (int * int) Collection.collection Collection.collection
