import copy
import math
import random

import pytest


class TSVReporterPlugin:
    def __init__(self):
        self.results = []
        self.processed_items = set()
        self.all_items = {}  # Stores test items by node ID

    def pytest_sessionstart(self, session):
        self.results = []
        self.processed_items = set()
        self.all_items.clear()

    def pytest_collection_modifyitems(self, items):
        # Store all test items by their nodeid
        for item in items:
            self.all_items[item.nodeid] = item

    def pytest_runtest_logreport(self, report):
        # Skip non-parameterized tests and avoid duplicates
        if report.nodeid in self.processed_items:
            return

        # Get the test item from our stored collection
        item = self.all_items.get(report.nodeid)
        if not item or not hasattr(item, 'callspec'):
            return

        # Process results only for call phase and setup skips/fails
        if report.when == 'call' or (report.when == 'setup' and report.outcome in ('skipped', 'failed')):
            param_value = next(iter(item.callspec.params.values())) if item.callspec.params else None
            self.results.append({
                'argument': param_value,
                'result': report.outcome
            })
            self.processed_items.add(report.nodeid)

    def pytest_sessionfinish(self, session, exitstatus):
        # Write TSV with sorted results (optional)
        with open('test_results.tsv', 'w', encoding='utf-8') as f:
            f.write('argument\tresult\n')
            for entry in sorted(self.results, key=lambda x: str(x['argument'])):
                arg = str(entry['argument']).replace('\t', '\\t').replace('\n', '\\n')
                f.write(f"{arg}\t{entry['result']}\n")


# def pytest_configure(config):
#     config.pluginmanager.register(TSVReporterPlugin())


def pytest_addoption(parser):
    parser.addoption("--subset-size", type=int, default=1_000_000, help="Random test subset size for parameterized tests")
    parser.addoption("--seed", type=int, default=None, help="Random seed for reproducibility")


def pytest_generate_tests(metafunc):
    subset_size = metafunc.config.getoption("--subset-size")
    if subset_size == 0:
        return
    seed = metafunc.config.getoption("--seed")
    random.seed(seed)  # Seed affects all parameterizations consistently
    old_markers = list(metafunc.definition.iter_markers("parametrize"))
    subset_size = math.ceil(math.pow(subset_size, 1/len(old_markers)))
    new_markers = []

    # Process all parameterized arguments in the test
    for marker in old_markers:
        param_args = marker.args[0]
        if "," in param_args:
            param_args = param_args.split(",")  # Split multiple arg names
        original_values = list(marker.args[1])
        
        if len(original_values) <= subset_size:
            new_markers.append((
                param_args,
                original_values,
                marker.kwargs.get("ids"),
                marker.kwargs.get("indirect", False)
            ))
            continue  # Skip sampling if not needed
        
        sampled = random.sample(original_values, subset_size)

        # Re-parameterize with sampled values
        new_markers.append((
            param_args,
            sampled,
            marker.kwargs.get("ids"),
            marker.kwargs.get("indirect", False)
        ))
    metafunc.definition.own_markers = []
    for args, values, ids, indirect in new_markers:
        print(args)
        metafunc.parametrize(
            args,
            values,
            ids=ids,
            indirect=indirect,
        )
