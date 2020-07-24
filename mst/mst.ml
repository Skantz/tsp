#load "str.cma";;
let query_random_node max_n node_lookup =
  let n : int = Random.int max_n             in
  let x = (Hashtbl.find node_lookup n)       in
  match x with
  | None -> let () = Printf.printf "%d\n" n in
    let line = read_line () in
    let line = Str.split (Str.regexp_string " ") line |> List.map int_of_string in
    Hashtbl.add node_lookup n line
  | Some x -> () (* Hashtbl.find node_lookup n  *) (* Hashtbl.find node_lookup x *)


let print_random_number max_n =
  let n : int = Random.int max_n in
  Printf.printf "%d\n" n

let n = read_int ()
let eps = read_float () 
let max_w = read_int ();;

let node_lookup = Hashtbl.create 1000;;

print_random_number 15
