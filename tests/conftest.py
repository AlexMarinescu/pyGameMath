import pytest


def pytest_configure(config):
    config.addinivalue_line('markers', 'defect(id): regression for a confirmed implementation or mathematical finding')
    config.addinivalue_line('markers', 'contract_question(id): proposed requirement awaiting API approval')


def pytest_collection_modifyitems(items):
    for item in items:
        for marker in item.iter_markers('defect'):
            item.add_marker(pytest.mark.xfail(strict=True, reason=marker.args[0]))
        for marker in item.iter_markers('contract_question'):
            item.add_marker(pytest.mark.xfail(strict=True, reason='Unapproved contract: '+marker.args[0]))


@pytest.fixture(autouse=True)
def audit_defect_metadata(request, record_property):
    marker = request.node.get_closest_marker('defect')
    if marker:
        record_property('defect_id', marker.args[0])
    question = request.node.get_closest_marker('contract_question')
    if question:
        record_property('contract_question_id', question.args[0])
