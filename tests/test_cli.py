import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class CLI(unittest.TestCase):
    def run_cli(self,*args):
        return subprocess.run([sys.executable,'-m','hgtlong',*map(str,args)],capture_output=True,text=True)

    def test_invalid_threshold(self):
        p=self.run_cli('call','--evidence','missing','--outdir','unused','--min-support',0)
        self.assertNotEqual(p.returncode,0)
        self.assertIn('must be positive',p.stderr)

    def test_failed_output_does_not_replace_previous_file(self):
        from hgtlong.cli import write_tsv
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'events.tsv';path.write_text('previous\n')
            def rows():
                yield {'event_id':'one'}
                raise ValueError('interrupted extraction')
            with self.assertRaisesRegex(ValueError,'interrupted'):
                write_tsv(path,rows())
            self.assertEqual(path.read_text(),'previous\n')
            self.assertEqual(len(list(Path(d).iterdir())),1)

    def test_empty_output_schema(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d); (p/'empty.tsv').write_text('read\ta0\ta1\tb0\tb1\n')
            result=self.run_cli('call','--evidence',p/'empty.tsv','--outdir',p/'out')
            self.assertEqual(result.returncode,0,result.stderr)
            fields=(p/'out/events.tsv').read_text().strip().split('\t')
            self.assertIn('n_a0',fields)
            self.assertIn('status',fields)
            self.assertFalse(json.loads((p/'out/run.json').read_text())['comparator_used_for_calling'])

    def test_reject_noncollated(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'panel.tsv').write_text('contig\tmag\na\tA\n')
            rows=[f'{name}\t0\ta\t1\t60\t1000M\t*\t0\t0\t*\t*\tNM:i:0\n' for name in ['r1','r2','r1']]
            (p/'reads.sam').write_text(''.join(rows))
            result=self.run_cli('extract','--sam',p/'reads.sam','--panel',p/'panel.tsv','--out',p/'evidence.tsv')
            self.assertNotEqual(result.returncode,0)
            self.assertIn('not query-name collated',result.stderr)

    def test_parallel_matches_serial(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'panel.tsv').write_text('contig\tmag\na\tA\nb\tB\n')
            with open(p/'reads.sam','w') as f:
                for i in range(300):
                    f.write(f'r{i}\t0\ta\t1001\t60\t2000M\t*\t0\t0\t*\t*\tNM:i:0\n')
                    f.write(f'r{i}\t256\tb\t5001\t0\t500S1000M500S\t*\t0\t0\t*\t*\tNM:i:0\n')
            for workers in [1,2]:
                result=self.run_cli('extract','--sam',p/'reads.sam','--panel',p/'panel.tsv',
                                    '--out',p/f'out{workers}.tsv','--workers',workers)
                self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual((p/'out1.tsv').read_bytes(),(p/'out2.tsv').read_bytes())
            self.assertEqual(json.loads((p/'out1.tsv.stats.json').read_text()),
                             json.loads((p/'out2.tsv.stats.json').read_text()))

    def test_comparison_does_not_hide_orientation_conflict(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/'events.tsv').write_text('event_id\tstatus\tctgA\tctgB\ta0\ta1\tb0\tb1\torientation\n'
                                       'E\tread_supported\ta\tb\t100\t200\t300\t400\t+\n')
            (p/'hdmi.tsv').write_text('HGT_ID\tONT_query_Reference_ID\tONT_subject_Reference_ID\tONT_query_Start_1based\tONT_query_End_1based\tONT_subject_Start_1based\tONT_subject_End_1based\tq_start\tq_end\ts_start\ts_end\n'
                                     'H\ta\tb\t101\t200\t301\t400\t101\t200\t400\t301\n')
            result=self.run_cli('compare','--events',p/'events.tsv','--hdmi',p/'hdmi.tsv','--out',p/'comparison.tsv')
            self.assertEqual(result.returncode,0,result.stderr)
            summary=json.loads((p/'comparison.tsv.summary.json').read_text())
            self.assertEqual(summary['cumulative_boundary_tolerance']['0']['hdmi_events'],0)
            self.assertEqual(summary['tiers']['orientation_conflict']['hdmi_events'],1)


if __name__=='__main__':unittest.main()
