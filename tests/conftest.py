import pytest


def pytest_configure(config):
    config.addinivalue_line('markers', 'defect(id): regression for a confirmed audit finding')


def pytest_collection_modifyitems(items):
    for item in items:
        for marker in item.iter_markers('defect'):
            item.add_marker(pytest.mark.xfail(strict=True, reason=marker.args[0]))


@pytest.fixture(autouse=True)
def audit_defect_metadata(request, record_property):
    marker = request.node.get_closest_marker('defect')
    if marker:
        record_property('defect_id', marker.args[0])
