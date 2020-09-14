## Cellranger debug record
Jan 19, 2020
Huitian (Yolanda) Diao

### 0. Running environment
- cellranger/3.1.0(default)
- Linux 64 (Scripps HPC)

### 1. Error message
*Reproducible error
```
2020-01-19 07:25:49 [runtime] (failed)          ID.Exp391_count.SC_RNA_COUNTER_CS.SC_RNA_COUNTER._BASIC_SC_RNA_COUNTER.MARK_DUPLICATES

[error] Pipestance failed. Error log at:
Exp391_count/SC_RNA_COUNTER_CS/SC_RNA_COUNTER/_BASIC_SC_RNA_COUNTER/MARK_DUPLICATES/fork0/chnk007-u6b032474ce/_errors

Log message:
thread 'main' panicked at 'Found 0 or >1 features for confidently mapped read/pair NS500124:192:HNJ53BGXC:4:13506:1000:1000': cr_stage/src/cmd_mark_dups.rs:285stack backtrace:
   0: martian::martian_main::{{closure}}::h91b2cb6e4ec208c2 (0x7f2d5488a444)
             at /mnt/home/jenkins/workspace/cellranger-3.1.0-2.2.5-build/sake/modules/cellranger-cs/3.1.0/lib/rust/.cargo/git/checkouts/martian-rust-35615836cc90309f/56bdd8d/src/lib.rs:418
   1: rust_panic_with_hook (0x7f2d5497c0a8)
             at src/libstd/panicking.rs:482
   2: std::panicking::begin_panic::h57e271f882aa496c (0x7f2d54826a69)
             at /rustc/fc50f328b0353b285421b8ff5d4100966387a997/src/libstd/panicking.rs:412
   3: cr_stage::cmd_mark_dups::process_barcode::hd751fc817f7edb3f (0x7f2d548348b9)
             at cr_stage/src/cmd_mark_dups.rs:285
   4: cr_stage::cmd_mark_dups::cmd_mark_dups::hd092a5adfde3f170 (0x7f2d5483a474)
             at cr_stage/src/cmd_mark_dups.rs:487
      <cr_stage::cmd_mark_dups::MarkDuplicatesStage as martian::MartianStage>::main::h974808f790ef4388
             at cr_stage/src/cmd_mark_dups.rs:571
   5: martian::do_main::h2e308a8e2c4f463b (0x7f2d54888de0)
             at /mnt/home/jenkins/workspace/cellranger-3.1.0-2.2.5-build/sake/modules/cellranger-cs/3.1.0/lib/rust/.cargo/git/checkouts/martian-rust-35615836cc90309f/56bdd8d/src/lib.rs:320
   6: martian::martian_main::h6fd98ef8e453a44e (0x7f2d5488a188)
             at /mnt/home/jenkins/workspace/cellranger-3.1.0-2.2.5-build/sake/modules/cellranger-cs/3.1.0/lib/rust/.cargo/git/checkouts/martian-rust-35615836cc90309f/56bdd8d/src/lib.rs:456
   7: cr_stage::main::hc37a84eb0c22f378 (0x7f2d54856ec6)
             at cr_stage/src/main.rs:84
   8: std::rt::lang_start::{{closure}}::h604fde077e847e29 (0x7f2d5482b3b2)
             at /rustc/fc50f328b0353b285421b8ff5d4100966387a997/src/libstd/rt.rs:64
   9: {{closure}} (0x7f2d5497b9c2)
             at src/libstd/rt.rs:49
      do_call<closure,i32>
             at src/libstd/panicking.rs:297
  10: __rust_maybe_catch_panic (0x7f2d5497f449)
             at src/libpanic_unwind/lib.rs:87
  11: try<i32,closure> (0x7f2d5497c4cc)
             at src/libstd/panicking.rs:276
      catch_unwind<closure,i32>
             at src/libstd/panic.rs:388
      lang_start_internal
             at src/libstd/rt.rs:48
  12: main (0x7f2d54857107)
  13: __libc_start_main (0x7f2d537f8d5c)
  14: <unknown> (0x7f2d548182c8)
  ```
  
### 2. Error message traceback
[Cellranger source codes](https://github.com/10XGenomics/cellranger)
1. "Found 0 or >1 features for confidently mapped read/pair" <br>
Match found: ellranger-master/lib/rust/cr_stage/src/cmd_mark_dups.rs
2. “grep "NS500124:192:HNJ53BGXC:4:13506:1000:1000" \*.fastq” <br>
**HTO** <br>
- Match found: HTO_S2_L004_I1_001.fastq, HTO_S2_L004_R1_001.fastq, HTO_S2_L004_R2_001.fastq <br>
- Match number: 44340 (14780 per file) <br>
- Total line count of HTO_S2_L004_R2_001.fastq: 32888052 (8222013 reads) <br>
- Match read percentage: 0.18% <br>
**cDNA** <br>
- Match found: cDNA_S2_L004_I1_001.fastq, cDNA_S2_L004_R1_001.fastq, cDNA_S2_L004_R2_001.fastq <br>
- Match number: 588651 (196217 per file) <br>
- Total line count of cDNA_S1_L004_I1_001.fastq: 436571308 (109142827 reads) <br>
- Match read percentage: 0.18% <br>

### 3. Strategy
1. Filter out "NS500124:192:HNJ53BGXC:4:13506:1000:1000" from fastq file & count: works for Lane 4
2. Fastqc input fastq files: fastq file quality is okay
3. Check bcl path: found Bcl file size incorrect. 
   - Error solved after re-uploading file and re-run `cellranger mkfastq` & `cellranger count`
   - No replicate ID in fastq files after re-run

### 4. Conclusion
1. Bcl file corrupted while uploading (internet problem: solved on 1/24/2020 by IT)
2. `cellranger mkfastq` suppress error message and created duplicated fastq IDs

## Future procautions:
1. Run check sum when uploading files
2. Check fastq file quality after `cellranger mkfastq`
