"""
Copyright 2026 Shun Fukagawa, Tsukasa Aso

  Licensed under the Apache License, Version 2.0.
  You may not use this file except in compliance with the License.
  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
"""

import sys
from dataclasses import dataclass
from typing import List, Dict

@dataclass
class DamageRecord:
    event_id: int
    dna_id: int
    chain_id: str
    residue_id: int
    compound_name: str
    atom_or_molecule_name: str
    damage_type: int  # 1: Direct, 2: Indirect

@dataclass
class SSBRecord:
    event_id: int
    chain_id: str
    residue_id: int
    damage_type: int  # 1: Direct, 2: Indirect

    def __eq__(self, other):
        if not isinstance(other, SSBRecord):
            return NotImplemented
        return (self.event_id == other.event_id and 
                self.chain_id == other.chain_id and 
                self.residue_id == other.residue_id)
    
    def __hash__(self):
        return hash((self.event_id, self.chain_id, self.residue_id))

class DamageAnalyzer:
    def __init__(self, input_file_name: str, base_pair_num: int):
        self.input_file_name = input_file_name
        self.base_pair_num = base_pair_num
        
        self.let_mean = 0.0
        self.let_std_dev = 0.0
        
        self.total_direct_damages = 0
        self.total_indirect_damages = 0
        self.total_damages = 0
        
        self.direct_ssbs = 0
        self.indirect_ssbs = 0
        self.total_ssbs = 0
        
        self.total_dsbs = 0
        
        self.indirect_to_total_ratio = 0.0
        self.ssb_to_dsb_ratio = 0.0
        
        self.damage_records: List[DamageRecord] = []
        self.ssb_records: List[SSBRecord] = []

    def get_input_file_name(self) -> str:
        return self.input_file_name

    def get_let_mean(self) -> float:
        return self.let_mean

    def get_let_std_dev(self) -> float:
        return self.let_std_dev
    
    def get_ssb_to_dsb_ratio(self) -> float:
        return self.ssb_to_dsb_ratio

    def analyze(self):
        # Get first row showing LET value
        self.parse_let_values(self.input_file_name)

        # Parse damage records
        self.parse_damage_records(self.input_file_name, self.damage_records)
        
        if not self.damage_records:
            self.no_result()
            return

        self.total_direct_damages = self.count_damages(self.damage_records, 1)
        self.total_indirect_damages = self.count_damages(self.damage_records, 2)
        self.total_damages = self.total_direct_damages + self.total_indirect_damages
        
        if self.total_damages > 0:
            self.indirect_to_total_ratio = self.total_indirect_damages / float(self.total_damages)
        else:
            self.indirect_to_total_ratio = 0.0
        
        # SSB calculation
        self.total_ssbs = self.count_ssbs(self.damage_records, self.ssb_records)
        
        # Direct and Indirect SSBs
        self.direct_ssbs = self.count_direct_ssbs(self.ssb_records)
        self.indirect_ssbs = self.count_indirect_ssbs(self.ssb_records)
        
        # DSB calculation
        self.total_dsbs = self.count_dsbs(self.ssb_records)
        
        if self.total_dsbs > 0:
            self.ssb_to_dsb_ratio = self.total_ssbs / float(self.total_dsbs)
        else:
            self.ssb_to_dsb_ratio = 0.0

    def parse_let_values(self, input_file_name: str):
        try:
            with open(input_file_name, 'r') as f:
                # Find the line with LET info
                for line in f:
                    if "LET =" in line:
                        self.let_mean = self.parse_let_mean(line)
                        self.let_std_dev = self.parse_let_std_dev(line)
                        break
        except Exception as e:
            print(f"Error parsing LET values: {e}", file=sys.stderr)

    def parse_let_mean(self, line: str) -> float:
        # Example line format: "# LET = 10.5 +- 0.5 keV / um"
        try:
            pos1 = line.find("LET = ")
            pos2 = line.find(" +- ")
            if pos1 == -1 or pos2 == -1:
                return 0.0 # Or raise error
            
            let_value_str = line[pos1 + 6 : pos2]
            return float(let_value_str)
        except ValueError:
            return 0.0

    def parse_let_std_dev(self, line: str) -> float:
        # Example line format: "# LET = 10.5 +- 0.5 keV / um"
        try:
            pos1 = line.find("+- ")
            pos2 = line.find(" keV")
            if pos1 == -1 or pos2 == -1:
                return 0.0
            
            let_std_dev_str = line[pos1 + 3 : pos2]
            return float(let_std_dev_str)
        except ValueError:
            return 0.0

    def parse_damage_records(self, input_file_name: str, records: List[DamageRecord]):
        try:
            with open(input_file_name, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split()
                    if len(parts) >= 7:
                        # eventID dnaID chainID residueId compoundName atomOrMoleculeName damageType
                        try:
                            record = DamageRecord(
                                event_id=int(parts[0]),
                                dna_id=int(parts[1]),
                                chain_id=parts[2],
                                residue_id=int(parts[3]),
                                compound_name=parts[4],
                                atom_or_molecule_name=parts[5],
                                damage_type=int(parts[6])
                            )
                            records.append(record)
                        except ValueError:
                            continue
        except Exception as e:
            print(f"Error reading records: {e}", file=sys.stderr)

    def count_damages(self, records: List[DamageRecord], damage_type: int) -> int:
        count = 0
        for record in records:
            if record.damage_type == damage_type:
                count += 1
        return count

    def count_ssbs(self, records: List[DamageRecord], ssb_records: List[SSBRecord]) -> int:
        for record in records:
            ssb = SSBRecord(record.event_id, record.chain_id, record.residue_id, record.damage_type)
            if ssb not in ssb_records:
                ssb_records.append(ssb)
        
        return len(ssb_records)

    def count_dsbs(self, ssb_records: List[SSBRecord], bp_distance: int = 10) -> int:
        dsb_count = 0
        
        # Group SSBs by eventID
        event_ssb_map: Dict[int, List[SSBRecord]] = {}
        for ssb in ssb_records:
            if ssb.event_id not in event_ssb_map:
                event_ssb_map[ssb.event_id] = []
            event_ssb_map[ssb.event_id].append(ssb)
            
        for event_id, ssbs in event_ssb_map.items():
            used = [False] * len(ssbs)
            
            for i in range(len(ssbs)):
                if used[i]:
                    continue
                
                for j in range(i + 1, len(ssbs)):
                    if used[j]:
                        continue
                    
                    ssb1 = ssbs[i]
                    ssb2 = ssbs[j]
                    
                    cond1 = (ssb1.chain_id != ssb2.chain_id)
                    
                    dist_val = abs(ssb1.residue_id - abs(ssb2.residue_id - self.base_pair_num))
                    cond2 = (dist_val <= bp_distance)
                    
                    if cond1 and cond2:
                        dsb_count += 1
                        used[i] = True
                        used[j] = True
                        break # Move to next SSB (i) loops
                        
        return dsb_count

    def count_direct_ssbs(self, ssb_records: List[SSBRecord]) -> int:
        count = 0
        for ssb in ssb_records:
            if ssb.damage_type == 1:
                count += 1
        return count

    def count_indirect_ssbs(self, ssb_records: List[SSBRecord]) -> int:
        count = 0
        for ssb in ssb_records:
            if ssb.damage_type == 2:
                count += 1
        return count

    def print_results(self):
        def print_line(label, value):
            print(f"{label:>35}   {value}")

        print_line("File Name: ", self.input_file_name)
        print_line("LET: ", f"{self.let_mean} +- {self.let_std_dev} keV/um")
        print_line("Total Direct Damages: ", self.total_direct_damages)
        print_line("Total Indirect Damages: ", self.total_indirect_damages)
        print_line("Total SSBs: ", self.total_ssbs)
        print_line("Direct SSBs: ", self.direct_ssbs)
        print_line("Indirect SSBs: ", self.indirect_ssbs)
        print_line("Total DSBs: ", self.total_dsbs)
        print_line("Indirect to Total Damage Ratio:(%) ", self.indirect_to_total_ratio * 100)
        print_line("SSB to DSB Ratio: ", self.ssb_to_dsb_ratio)
        
    def no_result(self):
        self.total_direct_damages = 0
        self.total_indirect_damages = 0
        self.direct_ssbs = 0
        self.indirect_ssbs = 0
        self.total_ssbs = 0
        self.total_dsbs = 0
        self.indirect_to_total_ratio = 0.0
        self.ssb_to_dsb_ratio = 0.0

def output_results_to_file(analyzers: List[DamageAnalyzer], output_file_name: str):
    try:
        with open(output_file_name, 'w') as f:
            for analyzer in analyzers:
                f.write(f"{analyzer.get_let_mean()} {analyzer.get_ssb_to_dsb_ratio()} {analyzer.get_let_std_dev()} 0\n")
    except Exception as e:
        print(f"Could not open output file {output_file_name}: {e}", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Please specify the input damage file (.txt).")
        print("Usage example: python3 print_damage_summary.py data/strand_breaks/StrandBreakInfo.root")
        exit()

    input_file_names = sys.argv[1:]

    try:
        base_pair_num_str = input("Enter the base pair number of the DNA to be analyzed: ")
        base_pair_num = int(base_pair_num_str)
    except ValueError:
        print("Invalid base pair number provided.")
        sys.exit(1)
    
    res = input("Do you want to output SSB/DSB ratio for each LET to a file? [y/n]: ")
    if res == "y":
        output_file_name = input("Enter the name of output file: ")
    
    analyzers = []
    
    for input_file_name in input_file_names:
        analyzer = DamageAnalyzer(input_file_name, base_pair_num)
        analyzer.analyze()
        analyzer.print_results()
        analyzers.append(analyzer)
        
    if output_file_name:
        output_results_to_file(analyzers, output_file_name)
